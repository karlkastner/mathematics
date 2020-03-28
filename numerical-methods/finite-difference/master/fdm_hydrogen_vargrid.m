% Wed Feb 29 19:33:22 MSK 2012
% Karl Kästner, Berlin
%
% [E N TA TE] = hydrogen_fdm_vargrid(d, L0, k)
%
% compute energy levels of a (confined) hydrogen atom
% based on the Schrödinger equation with Coulomb potential
% with the finite difference method on a variable grid
% simple setup - nuclues at the centre of the domain
%
% input
% d  : dimension 1, 2 or 3 (default 3)
% L0 : size of the confining box in atomimc units (default 40)
% k  : number of eigenvalues to be computed (default 10)
%
% output
% E : first k eigenvalues in columns
%     first column: approximation with lowest number of grid points
% N : number of grid points (correspond to columns of E)
% TA : set up time
% TE : time to compute the eigenvalues
function [E N TA TE] = fdm_hydrogen_vargrid(d, L0, k, N, f_)
	
	% check input arguments
	if (nargin < 1 || isempty(d) )
		d = 3;	
	end
	if (nargin < 2 || isempty(L0) )
		L0 = 40;	
	end
	if (nargin < 3 || isempty(k) )
		k = 10;	
	end
	if (nargin < 5 || isempty(f_) )
		f_ = 0;	
	end
	
	% prepare computation depending on dimension
	switch (d)
		case {1}
			f = @hydrogen_1d
			% number of grid points per dimension
			if (nargin < 4 || isempty(N))
				N = 2.^(6:20);
			end
			% ARPACK option, do not change
			eig_mode = 'SM';
			% manual shift for the eigenvalue
			s = -0.51;
			% parameter for grid setup
			e=1/1.2;
			afunc = @afunc_1d;
		case {2}
			f = @hydrogen_2d
			if (nargin < 4 || isempty(N))
				N = 2.^(1:9)
			end
			eig_mode = 'SM';
			s = -2.01;
			e=1/pi;
			afunc = @afunc_2d;
		case {3}
			f = @hydrogen_3d
			if (nargin < 4 || isempty(N))
				N = 2.^(2:6)
			end
			eig_mode = 'SA'
			s = -0.51;
			e=1/pi^2;
			afunc = @afunc_3d;
	end % switch d
	
	% preallocate matrices
	E = zeros(k,length(N));
	TA = zeros(1,length(N));
	TE = zeros(1,length(N));
	
	% compute eigenvalues for different grid sizes
	for idx=1:length(N)
		N
		n=N(idx);
		tic
		% set up the PDE
		[A I] = feval(f,n,L0,e);
		TA(idx)=toc()
		% enforce symmetrie (compensate for round off errors)
		A = 0.5*(A+A');
		% compute eigenvalues
			opts.isreal=1;
			opts.issym=1;
			opts.maxit = 3*n^d;
	%		opts.disp=2;
		switch (f_)
			case {0}
				eig_mode='SM'; % SA fails, SM uses shift and invert
				eig_mode='SA';
	%			s = 0;
				E(1:min(k,n^d),idx) = sort(real(eigs(A-s*I,[],min(k,n^d),eig_mode,opts)+s));
			case {1}
				B=speye(size(A));
				E(1:min(k,n^d),idx) = sort(eigs(@(x) minres(A-s*B,x,1e-12,n^d),n^d,min(k,n^d),s,opts));
			case {2}
				A = A - s*I;
				if (1 == d)
					A_ = full([[diag(A,-1);0] diag(A) [0;diag(A,+1)]]);
				else
					A_ = full([[diag(A,-n); zeros(n,1)] [diag(A,-1); 0] diag(A) [0; diag(A,+1)] [zeros(n,1); diag(A,+n)]]);
				end
				E(1:min(k,n^d),idx) = sort(eigs(@(x) afunc(A_,n,x),n^d,min(k,n^d),'SA',opts))+s;
		end
		% display eigenvalues
		E(1:k,:)
		TE(idx)=toc()
	end % for idx
	
end % function hydrogen_fdm_vargrid

% Grid and Laplacian setup in 1D
function [A I X L] = hydrogen_1d(n, L0,e)
	I = speye(n);
	h = L0/(n+1);
	X = h*(0:n+1) - L0/2;
	X = (X(end)/(e*(exp(e*X(end))-1))) * e * sign(X).*(exp(e*abs(X)) - 1);
%	c=1; p = 4;	% 3 for 1d, exp for 2d, 4 for 3d
%	X = 0.5*L0/X(end)^p*X.^p.*sign(X.^(p-1));
	% laplacian setup
	[L D1] = d_vargrid(X',2,2);
	% similarity transform to preserve symmetrie
	X = X(2:end-1);
%	L = sqrt(D1)*L*sqrt(D1);
	% coulomb potential
	R = diag(sparse(abs(X)));
	% Schrödinger equation
	A = -(0.5*L + inv(R));
end % hydrogen_1d

% Laplacian setup in 2D
function [A2 I2] = hydrogen_2d(n,L0,e)
	[A1 I1 X L1] = hydrogen_1d(n, L0,e);
	L2 = kron(L1,I1) + kron(I1,L1);
	R2 = sqrt(kron(diag(sparse(X.^2)),I1) + kron(I1,diag(sparse(X.^2))));
	A2 = -(0.5*L2 + inv(R2));
	I2= speye(n^2);
end % hydrogn_2d

% Laplacian setup in 3D
function [A3 I3] = hydrogen_3d(n,L0,e)
	[A1 I1 X L] = hydrogen_1d(n, L0,e);
	L3 = kron(L,kron(I1,I1)) + kron(kron(I1,L),I1) + kron(I1,kron(I1,L));
	R3 = sqrt(kron(kron(diag(sparse(X.^2)),I1),I1) ...
			+ kron(kron(I1,diag(sparse(X.^2))),I1) ...
			+ kron(I1,kron(I1,diag(sparse(X.^2)))));
	A3 = -(0.5*L3 + inv(R3));
	I3= speye(n^3);
end % hydrogen_3d

function y = afunc_1d(A,n,x)
	y = A(:,2).*x;
	y(2:end) = y(2:end) + A(1:end-1,1).*x(1:end-1);
	y(1:end-1) = y(1:end-1) + A(2:end,3).*x(2:end);
end % afunc_1d

function y = afunc_2d(A,n,x)
	y = A(:,3).*x;
	y(1:end-n+1) = y(1:end-n+1) + A(n:end,1).*x(n:end);
	y(2:end) = y(2:end) + A(1:end-1,2).*x(1:end-1);
	y(1:end-1) = y(1:end-1) + A(2:end,4).*x(2:end);
	y(n:end) = y(n:end) + A(1:end-n+1,5).*x(1:end-n+1);
end

function y = afunc_(A,n,x)
	y = A*x;
end


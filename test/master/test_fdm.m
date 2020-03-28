% Wed Feb 29 18:43:21 MSK 2012
% Karl Kästner, Berlin


function E = test_fdm(d)

	% number of computed eigenvalues
	k = 30;
	% dimension
	switch (d)
		case {1}
			N = 2.^(2:20);
			f = @hydrogen_1d;
			L0 = 40;
			k = 11;
			mode='SM' % SA
			s = -0.5;
			sigma = -0.5;
		case {2}
			N = 2.^(2:9);
			f = @hydrogen_2d;
			L0 = 80;
			s = -2;
			mode='SM'
			sigma = -2;
		case {3}
			N = 2.^(2:6);
			f = @hydrogen_3d;
			L0 = 40;
			mode='SA'
			s=0; %-1;
			sigma = -0.5;
	end % switch
			
	% allocate output arrays
	E = zeros(k,length(N));
	TE = zeros(1,length(N));
	TA = zeros(1,length(N));
	T = zeros(1,length(N));
	
	for idx=1:length(N)
		N
		n = N(idx)
	
		tic();
		[A I] = feval(f,n,L0);
		size(A)
		TA(idx) = toc()
		tic
		B = speye(n^d);
		opts.issym=1;
%		E(1:min(k,n^d),idx) = sort(eigs(@(x) minres(A-sigma*B,x,1e-7,1e3),n^d,min(k,n^d),sigma,opts))
%		E(1:min(k,n^d),idx) = sort(eigs(A,[],min(k,n^d),sigma,opts))
		E(1:min(k,n^d),idx) = sort(eigs(A-s*I,[],min(k,n^d),mode,opts))+s
		TE(idx) = toc
	end
	
	E_true = E(:,end);
	Err = E - E_true*ones(1,length(N))
	nErr = sum(Err.^2)
	loglog(N,nErr)
	
end % function test_fdm

function [A1 I1] = hydrogen_1d(n,L0)
	h = L0/(n+1);
	X = h*(1:n)' - L0/2;
	I1 = speye(n);
	L1 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	R1 = sqrt(diag(sparse(X.^2)));
	A1 = -(0.5*L1 + inv(R1));
end

function [A2 I2] = hydrogen_2d(n,L0)
	h = L0/(n+1);
	X = h*(1:n)' - L0/2;
	I = speye(n);
	L1 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	L2 = kron(L1,I) + kron(I,L1);
	R2 = sqrt(kron(diag(sparse(X.^2)),I) + kron(I,diag(sparse(X.^2))));
	A2 = -(0.5*L2 + inv(R2));
	I2= speye(n^2);
end

function [A3 I3] = hydrogen_3d(n,L0)
	h = L0/(n+1);
	X = h*(1:n)' - L0/2;
	I = speye(n);
	L1 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	L3 = kron(L1,kron(I,I)) + kron(I,kron(L1,I)) + kron(I,kron(I,L1));
	R3 = sqrt(kron(kron(diag(sparse(X.^2)),I),I) + kron(kron(I,diag(sparse(X.^2))),I) + kron(I,kron(I,diag(sparse(X.^2)))));
	A3 = -(0.5*L3 + inv(R3));
	I3= speye(n^3);
end

function y = afun(A,B,sigma,x)
	[y flag] = minres(A - sigma*B,x,1e-7,1e3);
end


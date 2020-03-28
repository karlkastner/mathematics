% Mon May 21 22:20:53 MSK 2012
% Karl Kästner, Berlin
%
% fdm_adaptive_grid(d, k, abstol, cnorm, L0, x0)
%
% compute energy levels of a (confined) hydrogen atom
% based on the Schrödinger equation with Coulomb potential
% with the finite difference method on a variable grid
% simple setup - nuclues at the centre of the domain
%
% input:
% d      : dimension 1, 2 or 3 (default 2)
% k      : number of eigenvalues to be computed (default 1)
% abstol : absolute tolerance (convergence criteria)
% cnorm  : norm to measure the error (1 == L1, 2=L3, otherwise L_oo)
%         note : only L1 norm converges         
% L0     : size of the confining box in atomimc units (default 10*ones(1,d))
% x0     : position of nucleus inside the cavity (default: 0.5*L)
%
% output:
% name : saves simulated data to a file identified with a time stamp and calls fdm_plot on the data
%
function name = fdm_adaptive_grid(d, k, abstol, cnorm, L0, x0)
	% ARPACK option, do not change
	eig_mode = 'SM';
	opts.isreal=1;
	opts.issym=1;
	n_max = 1e6/(2*d+1)*1; %/k;
	% there must be at least 5 points per axis to make this scheme work
	% and number of grid points must be even if atom in centre
	n0 = 6;

	% check input arguments
	if (nargin < 1 || isempty(d) )
		d = 2;	
	end
	if (nargin < 2 || isempty(k) )
		k = 1;	
	end
	if (nargin < 3 || isempty(abstol) )
		abstol = 1e-2;	
	end
	if (nargin < 4 || isempty(cnorm) )
		cnorm = 2;
	end
	% domain parameters
	if (nargin < 5 || isempty(L0) )
		L0 = 10*ones(1,d);
	end
	if (nargin < 6 || isempty(x0) )
		x0 = 0.5*L0;
	end
	E_true = [];

	% prepare computation depending on dimension
	switch (d)
		case {1}
			f = @hydrogen_1d;
			% manual shift for the eigenvalue
			s = -0.51;
			if (nargin() < 6)
				% multiplicity: 1
				E_true = -0.5./(1:10).^2;
				E_true = E_true(1:k)';
			end
		case {2}
			f = @fdm_schroedinger_2d;
			f_refine = @fdm_refine_2d;
			ref_mode = 1;
			% manual shift for the eigenvalue
			s = -2.01;
			% known eigenvalues for unconfined system
			if ( norm(x0 - 0.5*L0) < 1e-7 && min(L0) >= 10 )
				% multiplicity: 2*n - 1
				E_true = -0.5./([1 2*ones(1,3) 3*ones(1,5) 4*ones(1,7) 5*ones(1,9) 6*ones(1,11) 7*ones(1,13)]-1/2).^2;
				E_true = E_true(1:k)';
			end
		case {3}
			f = @fdm_schroedinger_3d;
			f_refine = @fdm_refine_3d;
			ref_mode = 2;
			% manual shift of the eigenvalue
			s = -0.51;
			if ( norm(x0 - 0.5*L0) < 1e-7 && min(L0) >= 10 )
				% multiplicity: n^2
				E_true = -0.5./[1 4*ones(1,4) 9*ones(1,9) 16*ones(1,16) 25*ones(1,25)];
				E_true = E_true(1:k)';
			end
		otherwise
			error('fdm_hydrogen_adaptive','here');
	end % switch d
	
	% preallocate matrices
	N = [];
	E = []; 
	Tr = [];
	err_est = [];

	% initial grid setup
	for idx=1:d
		X{idx} = linspace(0,L0(idx),n0)' - x0(idx);
	end
%		X{1} = linspace(0,L0(idx),6)' - x0(idx);
%		X{2} = linspace(0,L0(idx),8)' - x0(idx);
%		X{3} = linspace(0,L0(idx),10)' - x0(idx);

	kdx=1;
	idx=1;
	% while the number of found eigenvalues does not reach k
	while (kdx <= k)
		% set up the PDE
		tic();
		[A B D] = feval(f, X);
		Tr(1,idx) =toc();

		% record system size
		n = size(A,1);
		N(idx,1) = n;
		K(idx,1) = k; 
		disp(N)

		% compute eigenvalues
		opts.maxit = 3*n^d;
		k_ = min(k, n-1);
		tic();
		[v e] = eigs(A - s*B, [], k_, eig_mode, opts); % B is identity
		%[v e] = eigs(A, [], k_, s, opts); % B is identity
		Tr(2,idx) =toc();
		e = diag(e)+s;
		[e edx] = sort(e);
		v = v(:,edx);
		% undo similarity transform
		v = D*v; v=v/norm(v);
		E(1:k_,idx) = e;
		disp(E(1:k_,:));

		% mark and refine the mesh
		tic();
		% check how many eigenvalues have converged
		while (kdx <= k)
			[X_ err_max err] = feval(f_refine, X, v(:,kdx), ref_mode, cnorm);

			err_est(idx) = err_max;
			n_next = 1;
			for jdx=1:d
				n_next = n_next*size(X_{jdx},1);
			end
			% eigenvalue did not converge, refine the mesh
			if (err_max > abstol)
				if (n_next <= n_max)
					X = X_;
				end
				break;
			end
			% eigenvalue converged, proceed to the next eigenvalue
			disp(sprintf('eigenvalue %d converged',kdx))
			K(idx,1) = kdx;
			kdx=kdx+1;
		end
		Tr(3,idx) = toc();
		% display eigenvalues
		disp(err_est);
%		disp(E(1,:)-E_true(1));
		disp(Tr);
		if (n_next > n_max)
			disp(['Stopped as number of unknowns of next systems is ' num2str(n_next)]);
			break;
		end
		idx=idx+1;
	end % while kdx <= k

%	if (isempty(E_true))
%		E_true = E(:,end);
%	end

	% save data
	[void s] = system('date +%s'); s = regexprep(s, '\n', '');
	name = ['../dat-new/fdm-' s '.mat' ];
	disp(name);
	tag = 'fdm_vargrid';
	save(name, '-mat', 'tag', 'd', 'k', 'K', 'abstol', 'cnorm', 'L0', 'x0', 'N', 'Tr', 'E', 'err_est', 'E_true', 'v', 'err', 'X');
	fdm_plot(name);

end % function fdm_adaptive_grid()


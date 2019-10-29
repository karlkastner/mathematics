% Thu Feb 23 00:51:37 MSK 2012
% Karl Kästner, Berlin

function name = pdeeig_1d(afunc, bfunc, k, L0, x0, opt)

	%
	% check arguments
	%

	% Arpack options, do not change
	opts.isreal=1;
	opts.issym=1;
	eig_mode = 'SM';
	
	% assembly routine (1d)
	d = 1;

	if (nargin() < 1 || isempty(afunc))
		afunc = 1;
	end
	if (nargin() < 2)
		bfunc = [];
	end
	if (nargin() < 2 || isempty(k))
		k = 1;
	end
	if (nargin() < 3 || isempty(L0))
		L0 = 10;
	end
	if (nargin() < 4 || isempty(x0))
		x0 = 0;
	end
	if (nargin() < 5 || ~isfield(opt,'abstol'))
		opt.abstol = 1e-8;
	end
	if (nargin() < 5 || ~isfield(opt,'reltol'))
		opt.reltol = 1e-4;
	end
	if (nargin() < 5 || ~isfield(opt,'n_max'))
		opt.n_max = 1e6;
	end
	if (nargin() < 5 || ~isfield(opt,'order'))
		opt.order = 2;
	end
	if (nargin() < 5 || ~isfield(opt,'solver'))
		opt.solver = 'ARPACK';
	end
	if (nargin() < 5 || ~isfield(opt,'bcflag'))
		opt.bcflag=1;
	end
	if (nargin() < 5 || ~isfield(opt,'mflag'))
		opt.mflag=0;
	end
	if (nargin() < 5 || ~isfield(opt,'backward'))
		opt.backward=0;
	end
	if (nargin() < 5 || ~isfield(opt,'shift'))
		opt.shift=0;
	end
	if (nargin() < 5 || ~isfield(opt,'disp'))
		opt.disp = 0;
	end

	if (nargin() > 4 && isfield(opt,'E_true'))
		E_true = opt.E_true(:);
		E_true = E_true(1:k);
	end	

	int = @int_1d_gauss_4;

	%
	% initial mesh setup
	%
	tid = tic();
	n = 5;
	[P T Bc] = mesh_1d_uniform(n, L0);
	P = P - x0;


	%
	% main loop
	% 

	E = zeros(k,0);
	Tr = zeros(0,3);
	idx=1;
	kdx=1;
	% for each eigenvalue
	while (kdx <= k)
		% record number of triangle
		N(idx,2) = size(T,1);
		K(idx,1) = kdx;

		if (idx > 1) tid = tic(); end

		% todo, update neighbour in the refine routine
		Nm = neighbour_1d(P, T);

		% add additional triangle internal points for higher order methods
		%mesh = Mesh(P, T, Bc);
		% TODO promote 1D

		% record number of points
		N(idx,1) = size(P,1);

		% assemble the Laplacian matrix
		if (isscalar(afunc))
			A = afunc*assemble_1d_dphi_dphi(P, T, [], int);
		else
			A = assemble_1d_dphi_dphi(P, T, [], int);
		end

		% assemble the potential matrix
		if (isempty(bfunc))
			V = 0;
		elseif(isnumeric(bfunc))
			V = bfunc*assemble_1d_phi_phi(P, T, [], int);
		else
			V = assemble_1d_phi_phi(P, T, bfunc, int);
		end

		% combine to stiffness matrix
		A = A+V;

		% assemble the mass matrix
		B = assemble_1d_phi_phi(P, T, [], int);

		% apply boundary conditions
		[A B P__] = boundary_1d(A, B, Bc, opt.bcflag);
		Tr(idx,1) = toc(tid);
	
		% find the eigenvalues and eigenvectors
		k_ = min(k, size(P__,1));
		switch (opt.solver)
			case {'ARPACK_L'}
				% fix eigs bug
				if (1==size(A,1))
					v = 1;
					e=full(A/B);
				else
					% prepare the solver
					AsB = -(A - s*B);
					p = symamd(AsB);
					% todo, that is still not SPD
					L = chol(AsB(p,p));
					Lt = L';
					% solve
					[v e] = eigs(@(x) afunc(L, Lt, p, x), size(L,1), k_, 'SM', opts);
					e=(e+s);
				end
			case {'ARPACK'}
				% fix eigs bug
				if (1==size(A,1))
					v = 1;
					e=full(A/B);
				else
					[v e] = eigs(A-opt.shift*B,B,k_,eig_mode,opts);
					[e edx] = sort(diag(e));
					v = v(:,edx);
					e = e+opt.shift;
				end
			case {'JD'}
				try
					jdopts.maxiter=max(10,size(A,1));
		       	        	[v, ~, ~, ~, jds, jdt] = jdqz(A,B,k_,eig_mode,jdopts);
					e = diag(jds)/diag(jdt);
				catch cerr
					disp(cerr);
					dips('caught jd error')
					[v e] = eigs(A-opt.shift*B,B,k_,eig_mode,opts);
					e = e+opt.shift;
					pause
				end
				[e edx] = sort(diag(e));
				v = v(:,edx);
			case {'Lanczos'}
				e = lanczos_qr();
				e = sort(e');
				% eigenvectors are computed later
			otherwise
				st = dbstack();
				error(st.name,'Unknown Eigensolver');
		end
		E(1:k_,idx) = e;
		% do not update the shift - may converge to wrong eigenvalue
		%s = 1.1*e(1);
		Tr(idx,2) = toc(tid) - Tr(idx,1);
		
		% undo stripping of boundary points when boundary conditions where imposed strongly
		if (0 ~= opt.bcflag)
			v_ = zeros(size(P,1),k_);
			v_(P__,:) = v;
			v = v_;
		end


		% estimating the error and marking the cells for refinement
		% while kdx-eigenvalue has converged
		while (kdx <= k)
			if (1 == opt.backward) kdx_ = k-kdx+1; else kdx_ = kdx; end

			[h_side C] = regularity_1d(P,T,Bc);
			h_side_min(idx) = min(h_side);
			% estimate error of the next eigenvalue and mark cells for refinement
			[M err_est(idx,1) v_err] = mark_1d(P, T, v(:,kdx_), Nm, h_side, C);
			fprintf(1,'points: %d triangles %d value: %f error: %f\n', N(idx,1), N(idx,2), E(kdx_,idx), err_est(idx));

			% check for convergence
			if ( err_est(idx,1) > opt.abstol && err_est(idx,1) > opt.reltol*abs(E(kdx_,idx)) )
				break;
			end
			K(idx,1) = kdx;
			kdx=kdx+1;
		end
		% record number of elements to be refined
		MM(idx,1) = length(M);

		% stop if sufficiently many eigenvalues where found or grid became to large
		if (kdx > k || N(idx,1) > opt.n_max/2)
			Tr(idx,3) = toc(tid) - Tr(idx,1) - Tr(idx,2);
			fprintf(1,'assembly time: %f solver time: %f refinement time: %f\n',Tr(idx,1), Tr(idx,2), Tr(idx,3));
			break;
		end

		% refine the mesh
		[P T] = refine_1d(P, T, M);

		Tr(idx,3) = toc(tid) - Tr(idx,1) - Tr(idx,2);
		fprintf(1,'assembly time: %f solver time: %f refinement time: %f\n',Tr(idx,1), Tr(idx,2), Tr(idx,3));
	
		idx = idx+1;
	end % while k <= kdx

	if (~isfield(opt,'E_true'))
		E_true = E(:,end);
	end

	%
	% store data in a file
	%

	[void timestr] = system('date +%s'); timestr = regexprep(timestr, '\n', '');
%	name = ['../dat/fem-1d-' s '.mat'];
	name = [opt.folder '/fem-2d-' timestr '-' num2str(idx,'%02d') '.mat'];
	tag = 'fem_adaptive';
	save(name, '-mat', 'tag', 'd', 'k', 'L0', 'x0', ...
			'N', 'K', 'P', 'T', 'Bc', 'Nm', 'v', 'v_err', ...
			'E', 'E_true', 'err_est', 'Tr', 'MM', 'h_side_min', ...
			'opt');

	%
	% plot data
	%

	disp(name);
	if (opt.disp > 0)
		fem_plot_1d(name, 0);
	end

end % pdeeig_1d




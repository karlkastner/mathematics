% Fri Jul 27 15:43:01 MSK 2012
% Karl Kästner, Berlin
%
% function name = pdeeig_3d(afunc, bfunc, k, L0, x0, opt)
%
% find eigenvalues of PDE of the form:
%
%          a(x,y,z) \Delta u(x,y,z) + v(x,y,z) u(x,y,z) = \lambda u(x,y,z)
%
% afunc : coefficient of second derivative
%	  scalar or java function
% vfunc : coefficient of potential term
%
% k     : number of eigenvalues to compute
% L0    : domain size (box)
% x0    : origin of potential inside the domain
% opt   : structure with additional options
%  opt.E_true   : analyic eigenvalues (for test purposes)
%  opt.abstol   : absolute tolerance of computed eigenvalues
%  opt.backward : ignored in 3D at the moment
%  opt.bcflag   : controls implementation of dirichlet boundary condtions, 0 : weak, 1 : hard
%  opt.checkpoint : save intermediate results
%  opt.circular : ignored in 3D at the moment
%  opt.disp     : ignored at the moment
%  opt.folder   : output folder
%  opt.h_tol    : threshold for pregrading the grid at x0
%  opt.mflag    : 
%  opt.mode     : 0 : refine mesh such that lowest eigenvalue converges first and keep the refined mesh for the computation of subsequent eigenvalues
%                 1 : restart with a coarse mesh for each new eigenvalue
%  opt.n0       : number of grid points per axis in initial mesh
%  opt.n_max    : mximum number of mesh points
%  opt.poly     : degree of basis function polynomial
%  opt.reltol   : relative tolerance of computed eigenvalues
%  opt.shift    : shift of the eigenvalue solver
%  opt.solver   : eigenvalue solver
%  opt.t_max    : maximum run time
%
% output:
% name : saves simulated data to a file identified with a time stamp and calls fdm_plot on the data
%
function name = pdeeig_3d(afunc, vfunc, k, L0, x0, opt)

	%
	% check arguments
	%

	% Arpack options, do not change
	arpack_opt.isreal = 1;
	arpack_opt.issym = 1;
	eig_mode = 'SM';
	
	% dimension is 3
	d = 3;

	if (nargin() < 1 || isempty(afunc))
		afunc = 1;
	end
	if (nargin() < 2)
		vfunc = [];
	end
	if (nargin() < 2 || isempty(k))
		k = 1;
	end
	if (nargin() < 3 || isempty(L0))
		L0 = [10 10 10];
	end
	if (nargin() < 4 || isempty(x0))
		x0 = 0.5*L0;
	end
	if (nargin() < 5 || ~isfield(opt,'abstol'))
		opt.abstol = 1e-6;
	end
	if (nargin() < 5 || ~isfield(opt,'reltol'))
		opt.reltol = 1e-3;
	end
	if (nargin() < 5 || ~isfield(opt,'h_tol'))
		opt.h_tol = Inf;
	end
	if (nargin() < 5 || ~isfield(opt,'n_max'))
		opt.n_max = 8e4;
	end
	if (nargin() < 5 || ~isfield(opt,'poly'))
		opt.poly = 3;
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
	if (nargin() < 5 || ~isfield(opt,'n0'))
		opt.n0=[1 1 1]*3;
	end
	if (nargin() < 5 || ~isfield(opt,'backward'))
		opt.backward=0;
	end
	if (nargin() < 5 || ~isfield(opt,'shift'))
		opt.shift=0;
	end
	if (nargin() < 5 || ~isfield(opt,'circular'))
		opt.circular=0;
	end
	if (nargin() < 5 || ~isfield(opt,'folder'))
		opt.folder = '.';
		%opt.folder = '../dat-new/';
	end
	if (nargin() < 5 || ~isfield(opt,'t_max'))
		opt.t_max = Inf;
	end
	if (nargin() < 5 || ~isfield(opt,'disp'))
		opt.disp = 0;
	end
	if (nargin() < 5 || ~isfield(opt,'check'))
		opt.check = 0;
	end
	if (nargin() < 5 || ~isfield(opt,'checkpoint'))
		opt.checkpoint = 0;
	end
	if (nargin() < 5 || ~isfield(opt,'mode'))
		opt.mode = 0;
	end

	if (nargin() > 4 && isfield(opt,'E_true') && ~isempty(opt.E_true) && length(opt.E_true) >= k)
		E_true = opt.E_true(:);
		E_true = E_true(1:k);
	else
		E_true = [];
	end

	%
	% initial mesh setup
	%

	[void timestr] = system('date +%s'); timestr = regexprep(timestr, '\n', '');

	tid = tic();
	% create mesh-tree and generate the mesh
	[P T Bc] = mesh_3d_uniform(opt.n0, L0, x0);
	tree = Tree_3d(P, T, Bc);
	[v_sum0 a_sum0 h_eff_max volume area h_eff] = regularity_3d(P, T, Bc);
	convflag = 0;

	% points per element depending on basis functions
	% 4 10 20 35 ...
	pfac = nchoosek(opt.poly+d, d);

	if (isfield(opt,'int'))
		int = opt.int;
	else
	switch (opt.solver)
		case {'Lanczos'}
		switch (opt.poly)
			case {1} % linear basis functions
				int = @int_3d_nc_3;
			case {2} % quadratic basis functions
				int = @int_3d_nc_6;
			case {3} % cubic basis function
				int = @int_3d_nc_11;
			case {4} % quartic basis functions
				int = @int_3d_nc_20;
		end % switch opt.poly
		otherwise
		switch (opt.poly)
			case {1} % linear basis functions
				% do not choose the trapezoidal rule of origin is a grid point
				int = @int_3d_gauss_1;
			case {2} % quadratic basis functions         
				int = @int_3d_gauss_11;
			case {3} % cubic basis function
				int = @int_3d_gauss_24;
			case {4} % quartic basis functions 
				int = @int_3d_gauss_45;
		end % switch opt.poly
	end % opt.solver
	end

	%
	% variables
	%

	% computed eigenvalues, rows: sorted eigenvalues (ascending), coulmns: refinement iterations
	E = zeros(k,length(1));
	% measured run-time, rows: refinement iterations, columns: mesh-generation+promotion, assembly, eigensolve, refinement
	Tr = zeros(1,4);
	% target eigenvalue of each refinement iteration
	K = [];
	% number of mesh points in the mesh of each refinement iteration
	N = [];
	% computed eigenvectors of the last iteration
	v = [];
	% estimated error of the target eigenvector of the last iteration
	v_err = [];
	% estimated derivative norm of the targed eigenvector of the last iteration
	nH = [];
	% error norm of the target eigenvector of each refinement iteration
	err_est = [];
	% number of elements to be refined in each refinement iteration
	M_refine = [];
	% smallest element of the mesh of each refinement iteration
	h_min = [];
	% smallest measured value of the mesh degeneration of each refinement iteration
	degen = [];

	%
	% main loop
	%
	idx = 1;
	kdx = 1;
	k_arpack = min(2*k,k+5);
	% for each eigenvalue
	while (kdx <= k)
		K(idx,1) = kdx;

		if (idx > 1)
			tid = tic();
		end

		% generate the mesh from the tree
		mesh = tree.generate_mesh();
		mesh.element_neighbours();

%		if (opt.checkpoint)
%			save_data(idx);
%		end

		% promote the basis function to higher order (add points)
		if (opt.poly > 1)
			mesh.promote(opt.poly);
		end

		% record time for mesh generation and promotion
		Tr(idx,1) = toc(tid());

		mesh.prefetch();

		% record number of points
		N(idx,1) = mesh.np;
		N(idx,2) = mesh.nt;
		N(idx,3) = mesh.nb;

		% stop if number of unknowns exceeds the limit
		if (N(idx,1) > opt.n_max)
			idx = idx-1;
			break;
		end

		% assemble the Laplacian matrix
		if (isscalar(afunc))
			A = afunc*assemble_3d_dphi_dphi_java(mesh, [], int);
		else
			A = afunc*assemble_3d_dphi_dphi_java(mesh, afunc, int);
		end
		% assemble the potential matrix
		if (isempty(vfunc))
			V = 0;
		elseif(isnumeric(vfunc))
			V = vfunc*assemble_3d_phi_phi_java(mesh, [], int);
		else
			V = assemble_3d_phi_phi_java(mesh, vfunc, int);
		end
		% combine to stiffness matrix
		A = A+V;

		% assemble the mass matrix
		B = assemble_3d_phi_phi_java(mesh, [], int);

		% apply boundary conditions
		[P T Bc] = get_mesh_arrays(mesh);
		[A B P__] = boundary_3d(A, B, Bc, opt.bcflag);

		% record run time of assembly routines
		Tr(idx,2) = toc(tid);

		% find the eigenvalues and eigenvectors
		k_ = min(k, size(A,1));
		switch (opt.solver)
		case {'ARPACK'}
		% fix eigs bugs for "sparse" 1x1 matrices
			if (1 == size(A,1))
				v = 1;
				e = full(A/B);
			else
				flag = 1;
				while (0 ~= flag)
					arpack_opt.p = min(k_arpack,size(A,1));
					[v e flag] = eigs(A-opt.shift*B, B, k_, eig_mode, arpack_opt);
					if (0 ~= flag)
						k_arpack = k_arpack+5;
						fprintf(1,'increasing number of vectors to %d\n',k_arpack);
					end
				end
%			[v e] = eigs(A-opt.shift*B, B, k_, eig_mode, arpack_opt);
				[e edx] = sort(diag(e));
				v = v(:,edx);
				e = e+opt.shift;
			end
		case {'JD'}
			try
				jdopts.maxiter=max(10,size(A,1));
		       	       	[v e] = jdqz(A, B, k_, eig_mode, jdopts);
				%e = diag(jds)/diag(jdt);
			catch cerr
				disp(cerr);
				dips('caught jd error')
				[v e] = eigs(A-opt.shift*B,B,k_,eig_mode,arpack_opt);
				e = e+opt.shift;
			end
			[e edx] = sort(diag(e));
			v = v(:,edx);
			otherwise
				st = dbstack();
				error(st.name,'Unknown Eigensolver');
		end % switch
		k_ = length(e);
		E(1:k_,idx) = e;
		
	
		% undo stripping of boundary points when boundary conditions where imposed strongly
		if (0 ~= opt.bcflag)
			v__ = v;
			v_ = zeros(mesh.np, k_);
			v_(P__,:) = v;
			v = v_;
		end

		% record run time of eigensolver
		Tr(idx,3) = toc(tid);

		% prepare mesh operations	
		h_min(idx) = min(mesh.h_side);
		degen(idx) = min(mesh.degen);

		% estimating the error and marking the cells for refinement
		% while kdx-eigenvalue has converged
		while (kdx <= k)
			if (1 == opt.backward) kdx_ = k-kdx+1; else kdx_ = kdx; end
			kdx__ = min(max(1,kdx_),size(v,2));

			% calculate the pth-order partial derivatives
			dV  = mesh.dV(v(:, kdx__), opt.poly);
			% rate of convergence
			% rate = FEM.get_rate(1, opt.poly+1, 0);
			% rate = 2*opt.poly;
			rate = opt.poly+1;
			% calculate the error
			obj = mesh.estimate_error(dV, opt.poly+1, rate);
			% extract return values
			v_err = obj(1); err_est(idx,1) = obj(2); thresh = obj(3); nH = obj(4);

%{
			% error estimation by inner product
			[P T Bc] = get_mesh_arrays(mesh);
			vt = max( abs( [ v(T(:,1),kdx) v(T(:,2),kdx) v(T(:,3),kdx) ] ), [], 2);
			% h = double(mesh.h_max).^(opt.poly+1);
			h = double(mesh.h_max).^(2*opt.poly);
			err_est_(idx) = (h.*vt)'*nH/(v__(:,kdx)'*B*v__(:,kdx));
			v_err_ = (h.*abs(vt)).*nH;
			thresh_ = 0.5*max(v_err_);
%}

			% true error
			if (~isempty(E_true))
				err(idx) = abs(E(kdx,idx) - opt.E_true(kdx));
			else
				err(idx) = NaN;
			end

			% grade at the singularty if h_tol not yet reached
			if (h_min(idx) > opt.h_tol)
				pnorm = sum([sum(P(T(:,1),:).^2,2) sum(P(T(:,2),:).^2,2) sum(P(T(:,3),:).^2,2) sum(P(T(:,4),:).^2,2)],2);
				[vmin imin] = min(pnorm);
				M = imin';
			else
				% mark elements for refinement
				M = mark(v_err, thresh);
			end
			% display status
			fprintf(1,'number: %d points: %d triangles %d value: %f estimated error: %f\n', k, N(idx,1), N(idx,2), E(kdx,idx), err_est(idx));
			fprintf(1,'\tminimum edge %f mesh quality %f \n', h_min(idx), degen(idx));
%			fprintf(1,'\terror estimate 1: %f error estimate 2: %f true error: %f\n', err_est(idx), err_est_(idx), err(idx));

			% check for convergence or minimum number of eigenvalues
			% TODO use error estimate instead of the true error
			if ( kdx__ ~= kdx_ || (err(idx) > opt.abstol && err(idx) > opt.reltol*abs(E(kdx_,idx)) ))
			%if ( kdx__ ~= kdx_ || (err_est(idx,1) > opt.abstol && err_est(idx,1) > opt.reltol*abs(E(kdx_,idx)) ))
				break;
			end
			K(idx,1) = kdx;

			if (1 == opt.mode)
				% store the kth-eigenvalue
				E_a(kdx) = E(idx,kdx);
				% store the kth-eigenvector
				v_a(kdx).v = v(:,kdx);
				v_a(kdx).v_err = v_err;
				% store mesh for the k-th eigenvector
				mesh_s.P = mesh.P;
				mesh_s.T = mesh.T;
				mesh_s.Bc = mesh.Bc;
				mesh_s.N  = mesh.N;
				mesh_a(kdx) = mesh_s;
				% re-initialise the mesh
				[P T Bc] = mesh_3d_uniform(opt.n0, L0, x0);
				tree = Tree_3d(P, T, Bc);
			else
				mesh_a = [];
			end
			kdx=kdx+1;

		if (opt.checkpoint)
			save_data(idx);
		end
		end % while kdx <= k

		% record number of elements to be refined
		M_refine(idx,1) = length(M);

		% estimate number of points after refinement
		% for each refined tetra there are 7 = 8-1 additional tetras
		% each has pfac points, but those are mostly shared with each four neighbours
		n_next = round(N(idx,1) + size(M,1)*7*pfac/4);
		% stop if sufficiently many eigenvalues where found
		if (kdx > k)
			convflag = 1;
			break;
		end
 		% stop if grid became to large
		if (n_next > opt.n_max)
			fprintf(1,'Iteration was terminated before accuracy tolerance was met as the number of unknown exceeded the limit\n');
			break;
		end
		% stop if run time exceeded limit
		if (idx > 1 && sum(Tr(idx,1:3)+Tr(idx-1,4) > opt.t_max) )
			fprintf(1,'Iteration was terminated before accuracy tolerance was met as the run time of the last iteration exceeded the limit\n');
			break;
		end

		% refine the mesh
		switch (opt.mflag)
			case {0} % adaptive refinement
				tree.refine(M);
			case {1} % uniform refinement
				M=(1:size(T,1));
				tree.refine(M);
		end % switch opt.mflag

		% record run time of error estimation and refinement
		Tr(idx,4) = toc(tid);

		fprintf(1,'\tmesh generation time: %f assembly time: %f solver time: %f refinement time: %f\n', ...
				Tr(idx,1), Tr(idx,2)-Tr(idx,1), Tr(idx,3)-Tr(idx,2), Tr(idx,4)-Tr(idx,3));

		% verify the refined mesh
		if (opt.check)
			tree.check_neighbour();
			[P T Bc] = get_mesh_arrays(mesh);
			[v_sum a_sum h_eff_max volume area h_eff] = regularity_3d(P,T,Bc);
			da = (a_sum - a_sum0)/a_sum0;
			dv = (v_sum - v_sum0)/v_sum0;
			if ( abs(da) + abs(dv) > 1e-7)
				fprintf('Warning: Loss of accuracy volume: % e surface area %e \n', da, dv);
				%error('pdeeig_3d','here');
			end
		end % if opt.check

		if (opt.checkpoint)
			save_data(idx);
		end
		idx = idx+1;
	end % while k <= kdx

	% record run time of error estimation and refinement
	Tr(idx,4) = toc(tid);

	fprintf(1,'\tmesh generation time:  %f assembly time: %f solver time: %f refinement time: %f\n', ...
		               Tr(idx,1), Tr(idx,2)-Tr(idx,1), Tr(idx,3)-Tr(idx,2), Tr(idx,4)-Tr(idx,3));

	% run time of individual components
	Tr(:,2:end) = Tr(:,2:end) - Tr(:,1:end-1);

	if (isempty(E_true))
		E_true = E(:,end);
	end

	%
	% store data in a file
	%
	save_data(idx);

	%
	% plot data
	%

	disp(name);
	if (0 ~= opt.disp)
		fem_plot_3d(name, 0);
	end

	function save_data(idx)
		name = [opt.folder '/fem-3d-' timestr '-' num2str(idx,'%02d') '.mat'];
		tag = 'fem_adaptive'; % TODO, here uniform, pregraded, adatptive, h_tol
		% java-object to matlab-struct
		mesh_s.P  = mesh.P;
		mesh_s.T  = mesh.T;
		mesh_s.Bc = mesh.Bc;
		mesh_s.N  = mesh.N;

		save(name, '-mat', ...
				... % simulation parameters
				'tag', 'd', 'k', 'L0', 'x0', 'opt', ...
				... % iteration specific data
				'N', 'K', 'err_est', 'Tr', 'M_refine', 'h_min', 'degen', ...
				... % mesh-specific data
				'mesh_a', 'mesh_s', ...
				... % solution eigenvalue specific data
				'E_a', 'E', 'E_true', 'convflag', ...
				... % solution eigenvector specific data
				'v_a', 'v', 'v_err', 'nH' );
	end
end % pdeeig_3d()


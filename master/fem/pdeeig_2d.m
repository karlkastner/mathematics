% Mon May 21 11:30:15 MSK 2012
% Karl Kästner
%
% function name = pdeeig_2d(afunc, vfunc, k, L0, x0, opt)
%
% find eigenvalues of PDE of the form:
%
%          a(x,y) \Delta u(x,y) + b(x,y) u(x,y) = \lambda u(x,y)
%
% afunc : coefficient of second derivative
%	  scalar or java function
% vfunc : coefficient of potential term
%
% k     : number of eigenvalues to compute
% L0    : domain size (box)
% x0    : origin of potential inside the domain
% opt   : structure with additional options
%  opt.E_true     : analyic eigenvalues (for test purposes)
%  opt.abstol     : absolute tolerance of computed eigenvalues
%  opt.backward   : refine for highest eigenvalue first
%  opt.bcflag     : controls implementation of dirichlet boundary condtions, 0 : weak, 1 : hard
%  opt.checkpoint : save intermediate results
%  opt.circular   : project boundary points to the circle with radius L0
%  opt.disp       : 
%  opt.folder     : output folder
%  opt.h_tol      : threshold for pregrading the grid at x0
%  opt.mflag      : 
%  opt.n0         : number of grid points per axis in initial mesh
%  opt.n_max      : mximum number of mesh points
%  opt.poly       : degree of basis function polynomial
%  opt.reltol     : relative tolerance of computed eigenvalues
%  opt.shift      : shift of the eigenvalue solver
%  opt.solver     : eigenvalue solver
%  opt.t_max      : maximum run time
%
% output:
% name : saves simulated data to a file identified with a time stamp and calls fdm_plot on the data
%
function name = pdeeig_2d(afunc, vfunc, k, L0, x0, opt)

	%
	% check arguments
	%

	% Arpack options, do not change
	arpack_opt.isreal = 1;
	arpack_opt.issym = 1;
	eig_mode = 'SM';
	
	% dimension is 2
	d = 2;

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
		L0 = [10 10];
	end
	if (nargin() < 4 || isempty(x0))
		x0 = 0.5*L0;
	end
	if (nargin() < 5 || ~isfield(opt,'abstol'))
		opt.abstol = 1e-8;
	end
	if (nargin() < 5 || ~isfield(opt,'reltol'))
		opt.reltol = 1e-4;
	end
	if (nargin() < 5 || ~isfield(opt,'h_tol'))
		opt.h_tol = Inf;
	end
	if (nargin() < 5 || ~isfield(opt,'n_max'))
		opt.n_max = 6e5;
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
		opt.n0=[1 1]*3;
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
	% generate initial mesh
	[P T Bc X] = mesh_2d_uniform(opt.n0, L0, x0);
	mesh = Mesh_2d(P, T, Bc);
	[a_sum0 l_sum0 area l_bounday h_side s_angle C] = regularity_2d(P,T,Bc);
	convflag = 0;
	% inverse mapping : boundary -> triangle
	mesh.element_neighbours();
	Nm = mesh.N;

	if (isfield(opt,'int'))
		int = opt.int;
	else
	switch (opt.poly)
	case {1} % linear basis functions
		% (do not use the trapezoidal rule if gridpoint in origin)
		% int = @int_2d_nc_3; % trapezoidal rule
		int = @int_2d_gauss_3; % 1
	case {2} % quadratic basis functions
		int = @int_2d_gauss_6; % fourth order
	case {3} % cubic basis functions
		int = @int_2d_gauss_12; % sixth order
	case {4} % quartic basic functions
		int = @int_2d_gauss_16; % eigth order
	case {5} % quintic basic functions
		int = @int_2d_gauss_25; % 10th order
	otherwise
		st = dbstack();
		error(st.name,'Order of accuracy has to be an integer between 2 and 6');
	end % switch poly
	end

	%
	% main loop
	% 
	E = zeros(k,length(1));
	Tr = zeros(1,4);
	K = [];
	N = [];
	v = [];
	v_err = [];
	err_est = [];
	MM = [];
	h_min = [];
	nH = [];
	degen = [];
	idx = 1;
	kdx = 1;
	k_arpack = min(2*k,k+5);
	% for each eigenvalue
	while (kdx <= k)
		K(idx,1) = kdx;

		if (idx > 1)
			tid = tic();
		end

		if (opt.checkpoint)
			save_data(idx);
		end

		% add additional triangle internal points for higher order methods
		if (opt.circular)
			pdx = [Bc(:,1); Bc(:,2)];
			P(pdx,:) = project_circle(P(pdx,:),L0(1)/2,-(x0-L0/2));
		end
		mesh = Mesh_2d(P, T, Bc);
		mesh.N = Nm;

		% select the integration rule and refinement routine
		if (opt.poly > 1)
			mesh.promote(opt.poly);
		end

		% record time for mesh generation and promotion
		Tr(idx,1) = toc(tid);

		if (opt.circular)
			pdx = [mesh.Bc(:,1); mesh.Bc(:,2)];
			mesh.P(pdx,:) = project_circle(mesh.P(pdx,:),L0(1)/2,-(x0-L0/2));
		end

		% prefetcht the element matrices and precalculate element properties
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
			A = afunc*assemble_2d_dphi_dphi_java(mesh, [], int);
		else
			A = afunc*assemble_2d_dphi_dphi_java(mesh, afunc, int);
		end
		% assemble the potential matrix
		if (isempty(vfunc))
			V = 0;
		elseif(isnumeric(vfunc))
			V = vfunc*assemble_2d_phi_phi_java(mesh, [], int);
		else
			V = assemble_2d_phi_phi_java(mesh, vfunc, int);
		end
		% combine to stiffness matrix
		A = A+V;

		% assemble the mass matrix
		B = assemble_2d_phi_phi_java(mesh, [], int);

		% apply boundary conditions
		[P T Bc Nm] = get_mesh_arrays(mesh);
		[A B P__] = boundary_2d(A, B, Bc, opt.bcflag);

		% record run time of assembly routines
		Tr(idx,2) = toc(tid);
	
		% find the eigenvalues and eigenvectors
		k_ = min(k, size(A,1));
		switch (opt.solver)
			case {'ARPACK_L'}
				% fix eigs bug
				if (1==size(A,1))
					v = 1;
					e=full(A/B);
				else
					% prepare the solver
					AsB = -(A - opt.shift*B);
					p = symamd(AsB);
					% todo, that is still not SPD
					L = chol(AsB(p,p));
					Lt = L';
					% solve
					[v e] = eigs(@(x) afunc(L, Lt, p, x), size(L,1), k_, 'SM', arpack_opt);
					e=(e+opt.shift);
				end
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
			case {'Lanczos'}
				e = lanczos_qr();
				e = sort(e');
				% eigenvectors are computed later
			otherwise
				st = dbstack();
				error(st.name,'Unknown Eigensolver');
		end % switch
		E(1:k_,idx) = e;
		% do not update the shift - may converge to wrong eigenvalue
		%s = 1.1*e(1);
		
		% undo stripping of boundary points when boundary conditions where imposed strongly
		if (0 ~= opt.bcflag)
			v_ = zeros(size(P,1), k_);
			v_(P__,:) = v;
			v = v_;
		end

		% record run time of eigensolver
		Tr(idx,3) = toc(tid);


		% monitor mesh properties
		h_min(idx) = min(mesh.h_side);
		degen(idx) = min(mesh.degen);

		% estimating the error and marking the cells for refinement
		% while kdx-eigenvalue has converged
		while (kdx <= k)
			if (1==strcmpi('Lanczos',opt.solver))
				% do one shift and invert step to get the next eigenvector
				% TODO - use MINRES in higher dimensions
				v(:,kdx) = rand(size(v,1),1);
				v(:,kdx) = (A - e(kdx)*B) \ v(:,kdx);
			end
			if (1 == opt.backward)
				kdx_ = k-kdx+1;
			else
				kdx_ = kdx;
			end
			kdx__ = min(max(1,kdx_),size(v,2));

			% estimate (p+1)-th partial derivative
			dV  = mesh.dV(v(:,kdx__), opt.poly);
			% estimate the error
			obj = mesh.estimate_error(dV, 1, opt.poly+1, 0);
			v_err = obj(1); err_est(idx,1) = obj(2); thresh = obj(3); nH = obj(4);
			% mark elements for refinement
			M = mark(v_err, thresh);

			fprintf(1,'points: %d triangles %d value: %f estimated error: %f\n', N(idx,1), N(idx,2), E(kdx__,idx), err_est(idx));

			% check for convergence or minimum number of eigenvalues
			if ( kdx__ ~= kdx_ || (err_est(idx,1) > opt.abstol && err_est(idx,1) > opt.reltol*abs(E(kdx_,idx)) ))
				break;
			end
			K(idx,1) = kdx;
			kdx=kdx+1;
		end % while kdx <= k
		% record number of elements to be refined
		MM(idx,1) = length(M);


		% stop if sufficiently many eigenvalues where found
		if (kdx > k)
			convflag = 1;
			break;
		end
 		% stop if grid became to large
		if (N(idx,1) > opt.n_max)
			fprintf(1,'Iteration was terminated before accuracy tolerance was met as the number of unknown exceeded the limit\n');
			break;
		end
		% stop if run time exceeded limit
		if (idx > 1 && sum(Tr(idx,1:3))+Tr(idx-1,4) > opt.t_max )
			fprintf(1,'Iteration was terminated before accuracy tolerance was met as the run time of the last iteration exceeded the limit\n');
			break;
		end

		% stop if sufficiently many eigenvalues where found or grid became to large
%		if (kdx > k || N(idx,1) > opt.n_max || sum(Tr(max(1,idx-1),:)) > opt.t_max)
%			if (N(idx,1) > opt.n_max)
%				fprintf(1,'Iteration was terminated before accuracy tolerance was met as the number of unknown exceeded the limit');
%				convflag = 0;
%			elseif (sum(Tr(max(1,idx-1),:)) > opt.t_max)	
%				fprintf(1,'Iteration was terminated before accuracy tolerance was met as the run time of the last iteration exceeded the limit');
%				convflag = 0;
%			else
%				convflag = 1;
%			end
%			Tr(idx,3) = toc(tid) - Tr(idx,1) - Tr(idx,2);
%			fprintf(1,'assembly time: %f solver time: %f refinement time: %f\n',Tr(idx,1), Tr(idx,2), Tr(idx,3));
%			break;
%		end

		% refine the mesh
		switch (opt.mflag)
			case {0} % adaptive refinement
				[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, M);
				mesh = Mesh_2d(P,T,Bc);
				mesh.N = Nm;
			case {1} % uniform refinement
				M=(1:size(T,1));
				[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, M);
				mesh = Mesh_2d(P,T,Bc);
				mesh.N = Nm;
			otherwise % TODO: structured refinement with constant edge length
				error('pdeeig_2d','unknown refinement option');
		end % switch opt.mflag

		% record run time of error estimation and refinement
		Tr(idx,4) = toc(tid);

		fprintf(1,'\tmesh generation time: %f assembly time: %f solver time: %f refinement time: %f\n', ...
				Tr(idx,1), Tr(idx,2)-Tr(idx,1), Tr(idx,3)-Tr(idx,2), Tr(idx,4)-Tr(idx,3));
	
		% verify the refined mesh
		if (opt.check)
			% verify domain area and boundary length
			[P T Bc Nm] = get_mesh_arrays(mesh);
			[a_sum l_sum area l_boundary h_side s_angle C] = regularity_2d(P, T, Bc);
			da = (a_sum - a_sum0)/a_sum0;
			dl = (l_sum - l_sum0)/l_sum0;
			if ( abs(da) + abs(dl) > 1e-7)
				fprintf('Warning: Loss of accuracy volume: % e surface area %e \n', da, dl);
				%error('pdeeig_3d','here');
			end

			% verify element neighbours
			mesh = Mesh_2d(P, T, Bc);
			mesh.element_neighbours();
			Nm_ = double(mesh.N);
			if (norm(Nm - Nm_) > 0)
				fprintf('Warning: refinement routine did not properly update the neighbourhood relations\n');
				Nm = Nm_;
			end

		%if ( ~opt.circular && (abs(sum(0.5*mesh.determinant)-prod(L0)) > 1e-7 || abs(sum(mesh.l_boundary)-2*sum(L0)) > 1e-7))
		%	disp(sprintf('domain area: %f boundary length: %f\n', sum(0.5*mesh.determinant), sum(mesh.l_boundary)));
		%	st = dbstack();
		%	error(st.name,'inconsistent mesh');
		%end
		end % if opt.check

		idx = idx+1;
	end % while k <= kdx

	% record run time of last error estimation and refinement
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
		fem_plot_2d(name, 0);
	end

	function save_data(idx)
		[P T Bc Nm] = get_mesh_arrays(mesh);
		name = [opt.folder '/fem-2d-' timestr '-' num2str(idx,'%02d') '.mat'];
		tag = 'fem_adaptive';
		save(name, '-mat', 'tag', 'd', 'k', 'L0', 'x0', ...
				'N', 'K', 'P', 'T', 'Bc', 'Nm', 'v', 'v_err', ...
				'E', 'E_true', 'err_est', 'Tr', 'MM', 'h_min', 'nH', ...
				'degen', 'opt', 'convflag');
	end
end % function pdeeig_2d

function x = afunc(L,Lt,p,x)
	x(p) = L \(Lt\x(p));
end


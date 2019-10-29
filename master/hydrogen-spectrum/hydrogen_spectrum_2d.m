% Tue May 29 17:43:54 MSK 2012
% Karl Kästner, Berlin

function simulate_2d()
	format compact

	opt.t_max = 100;
	opt.t_max = 1000;
	opt.backward = 0;

	% generate a series of meshes for the solver run time comparison
	tic()
	k = 1;
	opt.poly = 1;
	L0 = 20*[1 1 1];
	opt.E_true   = -2;
	opt.abstol   = 1e-4; 
	opt.folder   = '../dat/2d-mesh-series';
	opt.reltol   = 1e-4;
	opt.shift    = -2.2;
%	opt.checkpoint = 1;
	mkdir(opt.folder);
	f = javaObject('Potential_2D_Coulomb');
	x0 = 0.5*L0 + opt.reltol*pi/3;
	opt.int = @int_2d_nc_3;
%	pdeeig_2d(-0.5, f,  k, L0, x0, opt);
	opt = rmfield(opt,'int')

	% compare convergence of different basis functions - one eigenvalue
	tic();
	k = 1;
	poly         = [1 2 3 4 5]';
	h_tol        = [];
	L0           = ones(length(poly),1)*[1 1]*128;
	x0           = 0.5*L0;
	opt.E_true   = -2;
	opt.abstol   = 1e-4; 
	opt.check    = 1;
	opt.circular = 0;
	opt.folder   = '../dat/2d-order-of-accuracy-1'
	opt.h_tol    = Inf;
	opt.reltol   = 1e-4;
	opt.shift    = -2.2;
	mkdir(opt.folder);
%	confinement_series_2d(k,L0,x0,poly,opt);
	T(1) = toc()

	% compare convergence of different basis functions - several eigenvalues
	opt.reltol   = 1e-3;
	poly = [1 2 3 4 5]';
	tic();
	k = 10;
	opt.E_true   = -0.5./([1 2*ones(1,3) 3*ones(1,5) 4*ones(1,7) ...
                                 5*ones(1,9) 6*ones(1,11) 7*ones(1,13)]-1/2).^2;
	opt.folder   = '../dat/2d-order-of-accuracy-10'
	mkdir(opt.folder);
%	confinement_series_2d(k,L0,x0,poly,opt);
	T(2) = toc()

% commom poperties for next runs
	k=16;
	h=0.25;
	opt.E_true = [];
	opt.abstol = 1e-4; 
	opt.poly   = 3;
	opt.reltol = 1e-3;

% square domain - unshifted
	tid = tic();
	L0  = 2.^(-1+h:h:7)'*[1 1];
	x0  = 0.5*L0;
	opt.circular   = 0;
	opt.folder = '../dat/2d-square-unshifted';
	mkdir(opt.folder);
	confinement_series_2d(k,L0,x0,[],opt);
	T(3) = toc(tid)

% square domain - shifted
	tid = tic();
	L0 = 32;
	x0=2.^(-1+h:h:7)'; x0 = 0.5*L0*[x0/max(x0) ones(size(x0))];
	L0 = L0*ones(size(x0));
	opt.circular   = 0;
	opt.folder = '../dat/2d-square-shifted';
	mkdir(opt.folder);
%	confinement_series_2d(k,L0,x0,[],opt);
	T(4) = toc(tid)

% circular domain - unshifted
	opt.check = 0;
	tid = tic();
	opt.abstol = 3.3e-4; % 1e-4
	L0  = 2.^(-1+h:h:7)'*[1 1];
	x0  = 0.5*L0;
	opt.circular = 1;
	opt.folder = '../dat/2d-circle-unshifted';
	mkdir(opt.folder);
%	confinement_series_2d(k,L0,x0,[],opt);
	T(5) = toc(tid)

% circular domain - shifted
	tid = tic();
	L0 = 32;
	x0=2.^(-1+h:h:7)'; x0 = 0.5*L0*[x0/max(x0) ones(size(x0))];
	%x0 = x0(11,:), pause
	L0 = L0*ones(size(x0));
	opt.circular = 1;
	opt.folder = '../dat/2d-circle-shifted';
	mkdir(opt.folder);
%	confinement_series_2d(k,L0,x0,[],opt);
	T(6) = toc(tid)
end

function Err = confinement_series_2d(k,L0,x0,poly,opt)
	f = Potential_2D_Coulomb;
	Err = [];
	for idx=1:length(L0)
		k
		fprintf('Domain size %f\n',L0(idx));
		fprintf('Position %f\n',x0(idx));
		if (~isempty(poly))
			poly(idx)
			opt.poly = poly(idx);
		end
		%obj = setfield(obj,'field',value)
%		try
			pdeeig_2d(-0.5, f,  k, L0(idx,:), x0(idx,:), opt);
%		catch err
%			Err(idx+22) = idx+22
%			err
%		end
	end
	Err
end


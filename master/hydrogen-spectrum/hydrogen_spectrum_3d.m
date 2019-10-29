% Sat Jul 28 11:27:12 MSK 2012
% Karl Kästner, Berlin

function simulate_3d()
	format compact

	opt.t_max = 100;
	opt.backward = 0;

	% generate a series of meshes for the solver run time comparison
	tic()
	k = 1;
	opt.poly = 1;
	L0 = 128*[1 1 1];
	opt.E_true   = -2;
	opt.abstol   = 1e-4; 
	opt.folder = '../dat/3d-mesh-series';
	opt.reltol   = 1e-4;
	opt.shift    = -2.2;
	opt.checkpoint = 1;
	mkdir(opt.folder);
	f = Potential_3D_Coulomb;
	x0 = 0.5*L0 + opt.reltol*pi/3;
	opt.int = @int_3d_nc_4;
	pdeeig_3d(-0.5, f,  k, L0, x0, opt);

return
	% compare convergence of different basis functions - one eigenvalue

	tic();
	k = 1;
	poly = [1 2 3 4]';
%	h_tol = [0.5 0.25 1/6 0.125];
%	h_tol = 0.5*[ 1 1 1 1 ];
	h_tol = [];
	L0           = ones(length(poly),1)*[1 1 1]*128;
	x0           = 0.5*L0;
	opt.E_true   = -0.5;
	opt.abstol   = 1e-4; 
	opt.circular = 0;
	opt.folder   = '../dat/3d-order-of-accuracy-1'
	opt.reltol   = 1e-4; 
	opt.shift    = -0.51;
	opt.h_tol    = double(Inf);
	opt.check    = 1;
	mkdir(opt.folder);
	confinement_series_3d(k,L0,x0,poly,h_tol,opt);
	T(1) = toc()

return

	% compare convergence of different basis functions - several eigenvalues
	opt.t_max = 1000;
	poly = [1 2 3 4];
	tic();
	k = 10
	opt.reltol   = 1e-3;
	opt.E_true = -0.5./([1 4*ones(1,4) 9*ones(1,9) 16*ones(1,16) 25*ones(1,25) 36*ones(1,36) 49*ones(1,49)]).^2;
	opt.folder = '../dat/3d-order-of-accuracy-10'
	mkdir(opt.folder);
	confinement_series_3d(k,L0,x0,poly,h_tol,opt);
	T(2) = toc()
return

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
	opt.folder = '../dat/square-unshifted';
	mkdir(opt.folder);
%	series_2d(k,L0,x0,[],opt);
	T(3) = toc(tid)

% square domain - shifted
	tid = tic();
	L0 = 32;
	x0=2.^(-1+h:h:7)'; x0 = 0.5*L0*[x0/max(x0) ones(size(x0))];
	L0 = L0*ones(size(x0));
	opt.circular   = 0;
	opt.folder = '../dat/square-shifted';
	mkdir(opt.folder);
%	series_2d(k,L0,x0,[],opt);
	T(4) = toc(tid)

% circular domain - unshifted
	tid = tic();
	L0  = 2.^(-1+h:h:7)'*[1 1];
	x0  = 0.5*L0;
	opt.circular = 1;
	opt.folder = '../dat/circle-unshifted';
	mkdir(opt.folder);
%	series_2d(k,L0,x0,[],opt);
	T(5) = toc(tid)

% circular domain - shifted
	tid = tic();
	L0 = 32;
	x0=2.^(-1+h:h:7)'; x0 = 0.5*L0*[x0/max(x0) ones(size(x0))];
	L0 = L0*ones(size(x0));
%x0 = x0(23:end,:); % 8,9,11,14,21
%L0 = L0(23:end,:);
	opt.circular = 1;
	opt.folder = '../dat/circle-shifted';
	mkdir(opt.folder);
%series_2d(k,L0,x0,[],opt);
	T(6) = toc(tid)
end

function Err = confinement_series_3d(k,L0,x0,poly,h_tol,opt)
	f = Potential_3D_Coulomb;
	Err = [];
	for idx=1:length(L0)
%		k
		fprintf('Position %f %f %f\n',x0(idx,:));
		fprintf('Domain size %f %f %f\n',L0(idx,:));
		if (~isempty(poly))
			poly(idx)
			opt.poly = poly(idx);
		end
		if (~isempty(h_tol))
			opt.h_tol = h_tol(idx);
		end
		%obj = setfield(obj,'field',value)
	%	try
			pdeeig_3d(-0.5, f,  k, L0(idx,:), x0(idx,:), opt);
	%	catch err
	%		Err(idx+22) = idx+22
	%		err
	%	end
	end
	Err
end


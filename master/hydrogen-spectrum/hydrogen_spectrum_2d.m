% Tue May 29 17:43:54 MSK 2012
% Karl Kästner, Berlin

function simulate_2d()
	format compact

	opt.t_max = 100;
	opt.backward = 0;

	% generate a series of meshes for the solver run time comparison
	tic()
	k = 1;
	opt.poly = 1;
	L0 = 20*[1 1 1];
	opt.E_true   = -2;
	opt.abstol   = 1e-4; 
	opt.folder = '../dat/2d-mesh-series';
	opt.reltol   = 1e-4;
	opt.shift    = -2.2;
	opt.checkpoint = 1;
	mkdir(opt.folder);
	f = Potential_2D_Coulomb;
	x0 = 0.5*L0 + opt.reltol*pi/3;
	opt.int = @int_2d_nc_3;
	pdeeig_2d(-0.5, f,  k, L0, x0, opt);

return

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
	series_2d(k,L0,x0,poly,opt);
	T(1) = toc()

	% compare convergence of different basis functions - several eigenvalues
	opt.t_max = 1000;
	opt.reltol   = 1e-3;
	poly = [3 4 5]';
	tic();
	k = 10;
	opt.E_true   = -0.5./([1 2*ones(1,3) 3*ones(1,5) 4*ones(1,7) 5*ones(1,9) 6*ones(1,11) 7*ones(1,13)]-1/2).^2;
	opt.folder   = '../dat/2d-order-of-accuracy-10'
	mkdir(opt.folder);
	series_2d(k,L0,x0,poly,opt);
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

function confinement_series_1d()
	afunc = -0.5;
	bfunc = @potential_coulomb;
	k = 25;
	L0 = 0.5*2.^(0:0.125:6);
	x0 = 0;
	opt.reltol = 1e-4;
	opt.abstol = 1e-10;
	opt.shift  = -0.6;
	for idx=1:1:size(L0,1)
		L0(idx)
		if (L0(idx) > 60)
			pdeeig_1d(afunc, bfunc, k, L0(idx), x0, opt);
		end
	end
end

function Err = series_2d(k,L0,x0,poly,opt)
	f = Potential_2D_Coulomb;
	Err = [];
	for idx=1:length(L0)
		k
		x0(idx,:)
		L0(idx,:)
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
%{
	for idx=1:5 %2:6
		opt.poly = idx
		pdeeig_2d(afunc, bfunc, k, L0, x0, opt);
		%pause(10)
	end
	name_c=regexp(ls(opt.folder,'-1'),'\n','split'); fem_plot_2d_series(name_c,opt.folder,0)
%}	

%	opt.folder='../dat/order-of-accuracy-16-b'
%	mkdir(opt.folder);
%	for idx=4:5 %2:6
%		opt.poly = idx
%		pdeeig_2d(afunc, bfunc, k, L0, x0, opt);
%		pause(10)
%	end
%	name_c=regexp(ls(opt.folder,'-1'),'\n','split'); fem_plot_2d_series(name_c,opt.folder,0)

%{

% single eigenvalue - iterate order of accuracy
pdeeig_2d(3,  1, 5e-3, 1e6, [161 161], [], 2, [], 1, 0);
pdeeig_2d(3,  1, 5e-3, 1e6, [161 161], [], 3, [], 1, 0); % Crash for 160x160
pdeeig_2d(3,  1, 5e-3, 1e6, [161 161], [], 4, [], 1, 0);
pdeeig_2d(3,  1, 5e-3, 1e6, [161 161], [], 5, [], 1, 0);
pdeeig_2d(3,  1, 5e-3, 1e6, [161 161], [], 6, [], 1, 0);

% many eigenvalues
pdeeig_2d(3, 23, 1e-3, 1e6, [160 160], [], 4);

% strongly confined
pdeeig_2d(3,  9, 1e-3, 1e6, [ 10  10], [], 4);

% mildly confined and shifted
pdeeig_2d(3,  9, 1e-3, 1e6, [ 80  80], [1 40], 4);

%}



	k = 20;
	opt.poly = 3;
	L0 = 20*[1 1];
	x0 = 0.5*L0;
%	opt.E_true   = -2;
	opt.abstol   = 1e-4; 
	opt.folder = '../dat/2012-12-2';
	opt.reltol   = 1e-4;
	opt.shift    = -2.2;
%	opt.checkpoint = 1;
	mkdir(opt.folder);
	f = javaObject('Potential_2D_Coulomb');
%	opt.int = @int_2d_nc_3;
	pdeeig_2d(-0.5, f,  k, L0, x0, opt);


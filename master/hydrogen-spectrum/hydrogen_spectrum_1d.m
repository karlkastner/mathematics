% 2012 Jul 28 11:29 (MSK)
% Karl Kästner, Berlin

function confinement_series_1d()
	afunc = -0.5;
	bfunc = @potential_coulomb;
	%k = 25;
	k = 5;
	h = 0.125;
	L0 = 2.^(-1+h:h:7)';
	x0 = 0;
	opt.reltol = 1e-4;
	opt.abstol = 1e-10;
	opt.shift  = -0.6;
	opt.folder = '../dat/1d-confinement';
	mkdir(opt.folder);
	for idx=1:1:size(L0,1)
		L0(idx)
		pdeeig_1d(afunc, bfunc, k, L0(idx), x0, opt);
	end
end % confinement_series_1d()


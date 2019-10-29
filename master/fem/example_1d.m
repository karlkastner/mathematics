% Tue Jun 19 16:12:20 MSK 2012
% Karl Kästner, Berlin

function example_1d(mode)
	switch(mode)
		case {1}
			%
			afunc = 1;
			bfunc = [];
			k = 1;
			L0 = 1;
			x0 = 0.0; % does not matter	
			opt.reltol = 1e-4;
			pdeeig_1d(afunc, bfunc, k, L0, x0, opt);
		case {2}
			% hydrogen
			afunc = -0.5;
			bfunc = @potential_coulomb;
			k = 5;
			L0 = 40;
			x0 = 0;
			opt.reltol = 1e-4;
			opt.shift = -0.6;
			pdeeig_1d(afunc, bfunc, k, L0, x0, opt)
		case {3}
			
	end
end % example_1d


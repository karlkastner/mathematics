% 2022-01-06 17:20:01.004124047 +0100
% note : modes S^p are not defined for p>0

function fail = test_spectral_density_brownian_phase_mode()

	syms f f0 s
	[fc,Sc] = spectral_density_brownian_phase_mode(f0,s)

	S = spectral_density_brownian_phase(fc,f0,s)
	
	% check value of Sc
	res1fun = matlabFunction(S-Sc);
	f0_ = 1.1;
	s_  = 0.3;
	res(1) = rms(res1fun(f0_,s_));
	
	% check locaction of maximum, dS/df = 0
	S      = spectral_density_brownian_phase(f,f0,s);
	dS_df  = diff(S,f);
	dS_df  = simplify(subs(dS_df,f,fc));
	resfun2 = matlabFunction(dS_df,'Vars',{'f0','s'});
	res(2) = rms(resfun2(f0_,s_));
	
	fail = rms(res) > 1e-4;
	
end


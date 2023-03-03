% 2022-01-06 16:40:44.611191816 +0100

function fail = test_spectral_debsity_brownian_phase_scale()

	syms f f0 s;
	S = spectral_density_brownian_phase(f,f0,s);
	S = matlabFunction(S)
	IS = quad(@(f) S(f,0.5,1),0,1e5)

	res = IS-1;
	fail = abs(res) > 1e-4;
end


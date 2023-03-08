

	Scy = 2;

	sy = spectral_density_brownian_phase_across_mode2par(Scy)

	Scy_ = spectral_density_brownian_phase_across(0,sy)

	fail = abs(Scy_ - Scy) > sqrt(eps)


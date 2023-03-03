% Thu 13 Oct 15:27:00 CEST 2022
function [f0,s] = spectral_density_brownian_phase_reg2par(reg,fc)
	Sc = reg./fc;
	[f0,s] = spectral_density_brownian_phase_mode2par(fc,Sc);
end


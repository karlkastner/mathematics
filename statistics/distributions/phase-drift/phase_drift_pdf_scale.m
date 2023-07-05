% 2022-01-06 16:55:34.460434615 +0100
% Karl Kästner, Berlin
%
%% normalization scale of the brownian phase spectral density
%
function I = spectral_density_brownian_phase_scale(f0,s)
	p  = (s.^2./pi^3);
	I = 1/(2*pi)*f0./p*pi^2;
end


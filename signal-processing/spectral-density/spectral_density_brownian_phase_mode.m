% Thu  6 Jan 16:00:43 CET 2022
% Karl Kästner, Berlin
%
%% mode (maximum) of the spectral density of the fourier series with brownian phase
%
function [fc,Sc] = spectral_density_brownian_phase_mode(f0,s)
	if (issym(s))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	p  = pi*s.^2;
	fc = f0.*sqrt(2*sqrt(1 + p.^2) - p.^2 - 1);
	Sc = p./(2*pi*f0.*(sqrt(1+p.^2) - 1));
end


% Wed  1 Dec 18:55:18 CET 2021
% Karl Kästner, Berlin
%
%% spectral density of a fourier series where the phase undergoes brownian motion
%% with standard deviation s per unit distance
%
%function [S,I] = spectral_density_brownian_phase(fx,f0,s,normalize)
function [S,I] = spectral_density_brownian_phase(fx,f0,s,normalize)
	if (nargin()<4)
		normalize = 0;
	end
	if (issym(s))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	p = pi_*s.^2;
	% this is already analytically normalized
	S = 4*p.*(p.^2 + fx.^2./f0.^2 + 1)./(2*pi_*f0.*(4*p.^2 + (p.^2 + fx.^2./f0.^2 - 1).^2));


	switch (normalize)
	case {1} % nuermically
		% TODO, use quad here
		I = spectral_density_area(fx,S);
		S = S./I;
	end
end % spectral_density_brownian_phase


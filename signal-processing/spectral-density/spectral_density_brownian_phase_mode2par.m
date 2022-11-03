% Fri  7 Jan 12:44:47 CET 2022
%
%% transform mode to parameters of the brownian phase spectral density
%
% function [f0, s] = spectral_density_brownian_phase_mode2par(fc,Sc)
function [f0, s] = spectral_density_brownian_phase_mode2par(fc,Sc)

	%r = root(z^3 - z^2 + z*(4*Sc^2*fc^2*pi^2 - 2) - 4*Sc^2*fc^2*pi^2, z, 1)
	rp = [1, -1, (4*Sc^2*fc^2*pi^2 - 2), -4*Sc^2*fc^2*pi^2];
	r = roots3(rp)
	% choose real root
	r = r(1);
	p2 = (  pi^2*((4*Sc^2*fc^2*pi^2 - 1)*r^2 ...
	      - 12*Sc^2*fc^2*pi^2 + 16*Sc^4*fc^4*pi^4 ...
              + 2*r)/(16*Sc^2*fc^2) );
	s = sqrt(pi/sqrt(p2));
	f0 = fc./sqrt(2*(pi.^2*s.^4 + 1).^(1/2) - s.^4*pi^2 - 1);
end % spectral_density_brownian_phase_mode2par


% Fri 22 Apr 16:00:51 CEST 2022
% mode (maximum) of the sd of the 2d bp
% inverse of the normalization constant
function Sc = spectral_density_bandpass_2d_scale(fc,p)
	error('not valid')
	Sc = (2^(2/3)-1)./(2*pi*fc.^2*4.^p) ...
	     .* (gamma(2*p + 1/3).*(3*p - 2)) ...
	     ./ (gamma(p + 1).*gamma(p + 1/3));
end


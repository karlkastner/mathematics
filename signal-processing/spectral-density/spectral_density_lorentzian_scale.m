% Mon 10 Jan 16:19:43 CET 2022
% Karl Kästner, Berlin
%
%% normalization scale of the lorentzian spectral density
%
function [IS] = spectral_density_lorentzian_scale(fc,p)
	IS = fc.*(atan(p) + pi/2)./p;
end


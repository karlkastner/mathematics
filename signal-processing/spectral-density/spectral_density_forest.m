% Fri 18 Nov 21:40:05 CET 2022
% spectral density of a random forest stand
function S = spectral_density_forest(fx,a)
	S = 2/pi*a./(a^2 + fx.^2);	
end

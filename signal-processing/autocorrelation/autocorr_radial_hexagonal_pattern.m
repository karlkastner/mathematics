% Sat 31 Dec 20:37:16 CET 2022
% Karl Kastner, Berlin
function a = autocorr_radial_hexagonal_pattern(x,lc)
	a = besselj(0,2*pi*x/lc);
end


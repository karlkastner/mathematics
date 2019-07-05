% Wed 16 May 12:05:52 CEST 2018
% Karl Kastner, Berlin
%
%% fourier series coefficients of a square pulse
%
function [an,bn] = fourier_coefficient_square_pulse(n,L,w)
	an = 2./(n*pi).*sin(n*pi*w/L);
	an(0 == n) = w/L;
	bn = zeros(size(an));
end

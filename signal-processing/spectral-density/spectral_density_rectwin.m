% Tue 17 May 15:44:53 CEST 2022
% function [S,T] = sd_rectwin(L,n,Lw)
function [S,T] = spectral_density_rectwin(L,n,Lw)
	f = fourier_axis(L,n);
	% transfer function
	T = sinc(f.*Lw);
	% density
	S = abs(T).^2;
end


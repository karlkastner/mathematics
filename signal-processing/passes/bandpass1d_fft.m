% Thu 25 Nov 12:34:20 CET 2021
% Karl Kastner, Berlin
%
%% filter input vector with a spatial (two-sided) bandpass in fourier space
%
function y = bandpass1d_fft(y,fc,p,dx)
	if (isvector(y))
		y = cvec(y);
	end
	n = size(y,1);
	x = dx*(0:n-1)';
	fx = fourier_axis(x);
	normalize = -1;
	S = spectral_density_bandpass_discrete(fx,fc,p,dx,normalize,'f');
	y = real(ifft(((S).^(p/2)).*fft(y)));
end

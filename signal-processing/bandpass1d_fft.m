% Thu 25 Nov 12:34:20 CET 2021
function y = bandpass1d_fft(y,fc,p,dx)
	n = size(y,1);
	x = dx*(0:n-1)';
	fx = fourier_axis(x);
	S = spectral_density_bandpass_discrete(fx,fc,p,dx,'f');
	sum(S)/x(end)
	y = real(ifft(((S).^(p/2)).*fft(y)));
end

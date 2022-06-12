% Wed 12 Jan 16:08:31 CET 2022
%
%% filter the input vector with a cosine-shaped highpass in frequency space
%
function [y_,S] = highpass1d_fft_cos(y,fc,L,p);
	if (isvector(y))
		y = cvec(y);
	end
	n  = size(y,1);
	f  = fourier_axis(L,n);
	S  = spectral_density_highpass_cos(f,fc,p);
	y_  = ifft(sqrt(S).*fft(y));
end


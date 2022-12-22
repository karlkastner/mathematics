%% Mon 19 Dec 17:03:02 CET 2022
%% smooth (filter) the  2D image z with a circular disk of radius nf
%% apply periodic boundary conditions
function [zbar, nf2, w, w_] = circfilt2(z, nf)
	[w,nf2] = circwin(size(z), nf);

	% Fourier transform of window
	fw = fft2(w);

	% by definition, the ft is real (but not positive)
	fw = real(fw);

	% apply mean filter by convolution theorem
	zbar = ifft2(fw.*fft2(z));

	if (isreal(z))
		zbar = real(zbar);
	end
end


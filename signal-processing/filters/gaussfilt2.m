% Wed 21 Dec 11:27:03 CET 2022
%% smooth (filter) the  2D image z with a gaussian window
%% apply periodic boundary conditions
function [zbar, dof, w] = gaussfilt2(z, nf)
	
	[w, dof] = gausswin2(size(z), nf);

	% Fourier transform of window
	fw = fft2(w);

	% by definition, the ft is real
	fw = real(fw);

	% apply mean filter by convolution theorem
	zbar = ifft2(fw.*fft2(z));

	if (isreal(z))
		zbar = real(zbar);
	end
end


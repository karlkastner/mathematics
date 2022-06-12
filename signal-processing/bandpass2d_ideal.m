% Wed 27 Apr 10:55:29 CEST 2022
function [y,S,R,r]=bandpass2d_ideal(x,L,dx,varargin)
%	y = lowpass2d_ideal(x,L,varargin{:});
%	y = highpass2d_ideal(y,L,varargin{:});
	[S,R,r] = spectral_density_bandpass2d_ideal(size(x),dx,L,varargin{:});
	y = ifft2(sqrt(S).*fft2(x));
end

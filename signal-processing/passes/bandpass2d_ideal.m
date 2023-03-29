% Wed 27 Apr 10:55:29 CEST 2022
function [y,S,R,r]=bandpass2d_ideal(x,L,dx,varargin)
	[S,R,r] = bandpass2d_pdf_discrete(size(x),dx,L,varargin{:});
	y = ifft2(sqrt(S).*fft2(x));
end

% Fri 22 Apr 16:00:51 CEST 2022
% scale (normalization) factor of the spectral density of the 2D bandpass filter
% the scale factor is equal to the maximum of the density
% mode (maximum) of the sd of the 2d bp
% inverse of the normalization constant
function [Sc] = bandpass2d_pdf_mode(L,n,a,order)
	fr = fourier_axis(L,n);
	fr = fr(fr>=0);
	S  = bandpass2d_pdf(fr,a,order);
	df  = 1/L;
	ciS = sum(S)*df;
	Sc  = 1./ciS;
end


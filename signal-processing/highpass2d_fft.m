% Tue 20 Jul 13:53:52 CEST 2021
function [y] = highpass2d_fft(x,rho,varargin)
	y = x - lowpass2d_fft(x,rho,varargin{:});
end

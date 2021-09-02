% Sat 26 Jun 10:59:41 CEST 2021
function [y] = highpass1d_implicit(x,rho,varargin)
	y = x - lowpass1d_implicit(x,rho,varargin{:});
end

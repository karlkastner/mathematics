% Thu 24 Jun 16:27:03 CEST 2021
function [y] = highpass2d_implicit(x,rho,varargin)
	y = x - lowpass2d_implicit(x,rho,varargin{:});
end

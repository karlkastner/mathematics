% Sat 26 Jun 10:59:41 CEST 2021
function [y] = highpass1d_implicit(x,rho,order,varargin)
	% y = x - lowpass1d_implicit(x,rho,varargin{:});
	if (nargin()<3)
		order = 1;
	end
	y = x;
	for idx=1:order
		y = y - lowpass1d_implicit(y,rho,1,varargin{:});
	end
end


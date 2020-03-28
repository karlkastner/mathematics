% 2014-11-22 13:33:33.590366360 +0100
%TODO : at the lower end the cdf is actually known exact,
%       so why interpolating it to the centres?
function [H, obj] = cdf(obj)
	H = obj.cdfS(obj.h);
end


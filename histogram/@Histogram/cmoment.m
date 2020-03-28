% Di 8. Mär 10:05:07 CET 2016
% Karl Kastner, Berlin
% central moments
function [m obj] = cmoment(obj,order,varargin)
	m = cmomentS(obj.h,obj.edge,order,varargin{:});
end


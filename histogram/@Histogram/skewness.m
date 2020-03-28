% Mon Nov 24 11:41:30 CET 2014
% Karl Kastner, Berlin

function [sk obj] = skewness(obj)
	sk = obj.skewnessS(obj.h,obj.edge,obj.method);
end


% Do 18. Sep 14:15:19 CEST 2014
% Karl Kastner, Berlon

function [q obj] = quantile(obj,p)
	q = obj.quantileS(obj.h,obj.edge,p);
end % quantile


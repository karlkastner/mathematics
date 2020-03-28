% Mon Nov 24 11:47:01 CET 2014
% Karl Kastner, Berlin

function [kurt obj] = kurtosis(obj)
	kurt = obj.kurtosisS(obj.h,obj.edge,obj.method);
end % kurtosis


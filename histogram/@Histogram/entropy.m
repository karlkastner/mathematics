% So 21. Sep 11:04:35 CEST 2014
% Karl Kastner, Berlin

function [Esum Ebin obj] = entropy(obj)
	[Esum Ebin] = obj.entropyS(obj.h,obj.edge);
end


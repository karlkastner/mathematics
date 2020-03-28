% Do 18. Sep 12:08:05 CEST 2014
% Karl Kastner, Berlin

% mode (peak) of a histogram 
% H      : [nhist,nbin] set of historgram computed on identical grid
% centre : [1xnbin] histogram bin centres
function [mo obj] = mode(obj)
	mo = obj.modeS(obj.h,obj.edge);
end % histmode


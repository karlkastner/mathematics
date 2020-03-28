% Do 18. Sep 12:05:25 CEST 2014
% Karl Kastner, Berlin

% compute mean of a histogram
% H      : [nhist,nbin] set of historgram computed on identical grid
% centre : [1xnbin] histogram bin centres
function [mu obj] = mean(obj)
	mu = obj.meanS(obj.h,obj.edge,obj.method);
%	 n  = size(h,1);
%	 mu = sum(bsxfun(@times,h,rvec(centres)),2);
end % histmean


% 2015-07-31 17:51:12.505066364 +0200
% Karl Kastner, Berlin
%
%% class for quasi-distributions from a set of sampling points
classdef PDF < handle
	properties
		n;
		CDF
		rank
		val
		s
		m
	end
	methods
	function obj = PDF(x)
		fdx = isfinite(x);
		obj.n = sum(fdx);
		[obj.val, obj.rank] = sort(x(fdx));
		obj.CDF = (1:obj.n)/(obj.n+1);
		% make values unique
		[obj.val udx] = unique(obj.val);
		obj.CDF = obj.CDF(udx);
		q = quantile(x,normcdf([-1 0 1]));
		obj.m = q(2);
		obj.s = 0.5*(q(3)-q(1));
	end
	function cdf = cdf(obj,x)
		cdf = interp1(obj.val,obj.CDF,x,'linear');
		% limit
		cdf = min(max(cdf,0),1);
	end
	function x = normalise(obj,x)
		 x = obj.s*norminv(obj.cdf(x)) + obj.m;
	end
	end %methods
end % class PDF


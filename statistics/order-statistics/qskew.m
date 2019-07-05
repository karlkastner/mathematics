% 2015-08-03 14:23:32.800719870 +0200
%
%% skewness estimated from quantiles
%%
%% Note : this is a measurement of shape-symmetry and yields the same value for the
%%        skew-normal distribution as "skewness"
%%        However, this is an own statistic and hence requires different
%%        methods for calculating P-values and hypothesis testing
%	TODO, normalise, such that value and derivative at 0 are identical to skewness
%	for skew-normal distribution
function sk = qskew(x,scaleflag,p)
	if (nargin < 3)
		p = normcdf(-1);
	end
	p = min(p,1-p);
	q = quantile(x,[p 0.5 1-p]);
	sk = qskewq(q,p);
end


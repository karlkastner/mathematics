% Mon Oct 14 09:01:25 UTC 2013
% Karl Kästner, Berlin
%
%% Hotelling's T-squared cumulative distribution
%
function alpha = t2cdf(t2,p,n)
	% transform t-squared into F-distribution
	% f = t2*(nx + ny - p - 1)/((nx + ny - 2)*p);
	f = t2*(n - p + 1)/(n*p);
	alpha = fcdf(f, p, n - p + 1);
end % t2cdf()


% 2013-10-13 13:15:50 UTC
% Karl Kästner, Berlin
%
%% inverse of Hotelling's T-squared cumulative distribution
%
function t2 = t2inv(alpha, p, n)
	% get the f-distribution
	f  = finv(alpha, p, n - p + 1);
	% transform the f-distribution into the t^2-distribution
	t2 = f*(n*p)/(n-p+1);
end % t2inv()


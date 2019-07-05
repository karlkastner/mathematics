% 2014-12-11 18:23:39.088410624 +0100
% Karl Kastner, Berlin
%
%% weighted median regression 
%% c.f. Scholz, 1978
%
function [c, res2, res] = weighted_median_regression(x,y)
	n = length(x);
	[ x sdx] = sort(x);
	y = y(sdx);
	% pairwise slopes
	X  = repmat(x,1,n);
	dx = X - X';
	Y  = repmat(y,1,n);
	dy = Y - Y';
	b = dy./dx;
	mask = triu(true(n,n),+1);
	beta = median(b(mask));
	alpha = median(y - beta*x);
	c = [alpha; beta];
	res = y - alpha - beta*x;
	res2 = res'*res/(n-2);

	
%	dx = diff(dx);
%	dy = diff(y);
%	b = dy./
end


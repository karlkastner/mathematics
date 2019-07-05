% Thu Dec 11 17:56:38 CET 2014
% Karl Kastner, Berlin
%
%% gini regression
%
function [c, res2, res] = ginireg(x,y)
	[x sdx] = sort(x);	
	y = y(sdx);

	n    = length(x);
	idx  = (1:n-1)';
	dx   = diff(x);
	dy   = diff(y);
	w    = (n-idx).*idx.*dx;
	w    = w/sum(w);
	% slope
	beta  = w'*(dy./dx);
	% intercept
	alpha = mean(y) - beta*mean(x);
	% residual
	res = y - alpha - beta*x;
	res2 = (res'*res)/(n-2);
	c = [alpha; beta];
end % ginireg


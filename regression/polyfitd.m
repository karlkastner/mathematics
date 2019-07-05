% Thu  1 Sep 18:37:05 CEST 2016
% Karl Kastner, Berlin
%
%% fit a polynomial of order n to a set of sampled values and sampled values
%% of the derivative
%%
%% x0 must contain at least for conditioning as otherwise the intercept
%% cannot be determined
%
function [p,A] = polyfitd(x0,y0,x1,y1,n)
	x0  = cvec(x0);
	y0  = cvec(y0);
	xd1 = cvec(x1);
	A0  = vander_1d(x0,n);
	A1  = vanderd_1d(x1,n,1);
	A = [A0;A1];
	y = [y0;y1];
	p = A \ y;
	p = flipud(p);
end


% 2016-08-22 12:38:58.598026551 +0200
% Karl Kastner, Berlin
%% extrema of a polynomial
function [x y minflag realflag] = poly_extrema(p)
	% first derivative
	p1    = polyder(p);
	% zeros of the first derivative, x values of extrema and sadlle points
	x     = roots(p1);
	% second derivative
	p2        = polyder(p1);
	% values of the second derivative
	y2        = polyval(p2,x);
	% minima
	minflag   = (y2 > 0);
	% real roots
	realflag = abs(imag(x)) < sqrt(eps)*abs(x);
	% function value
	y = polyval(p,x);
end


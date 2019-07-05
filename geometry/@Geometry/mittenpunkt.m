% Sat 10 Sep 13:25:25 CEST 2016
%
%% mittenpunkt of a triangle
function [x0 y0] = mittenpunkt(x,y)
	xl = left(x);
	xr = right(x);
	yl = left(y);
	yr = right(y);

	l = hypot(xr-xl,yr-yl);
	ll = left(l);
	lr = right(l);

	p = l.*(ll+lr-l)

%	p = x.*(xl+xr-x)
%	q = y.*(yl+yr-y)

	x0 = p*x.';
	y0 = p*y.';
	
end

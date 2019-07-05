% 2016-09-10 13:34:41.087849833 +0200

%% nagelpoint of a triangle

function [x0 y0] = nagelpoint(x,y)
	xl = left(x);
	xr = right(x);
	yl = left(y);
	yr = right(y);
	l = hypot(xr-xl,yr-yl);
	s = 0.5*sum(l);
	p = s-l;
	x0 = p*x.';
	y0 = p*y.';
end


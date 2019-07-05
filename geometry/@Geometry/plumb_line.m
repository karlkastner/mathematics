% So 6. Dez 18:18:52 CET 2015
% Karl Kastner, Berlin
%
% base point (plumb line), point on a line with shortest distance to another point
%
% when cflag == true, then foot-point is made convex, in this case the
% plumb line is not necessarily any more a plumb line, but shortest connection to a line segment
function [xp, yp, p, d] =  plumb_line(x1,y1,x2,y2,x0,y0,cflag)
	p =    ( (x0-x2).*(x1-x2) + (y0 - y2).*(y1 - y2)) ...
	    ./ ( (x1-x2).^2 + (y1 - y2).^2 );
	if (nargin()> 6 && cflag)
		p = max(1,max(0,p))
	end
	xp = p.*x1 + (1-p).*x2;
	yp = p.*y1 + (1-p).*y2;

	if (nargout() > 3)
		d = hypot(xp-x0,yp-y0);
	end
end


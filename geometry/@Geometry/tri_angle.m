% 2016-09-25 15:02:43.699837127 +0200
%% cos of angles of a triangle
% function cosa = tri_angle(x,y)
function cosa = tri_angle(x,y)
	xc = x; %squeeze(xy(1,:,:));
	xl = left(xc);
	xr = right(xc);
	yc = y; %squeeze(xy(2,:,:));
	yl = left(yc);
	yr = right(yc);

	dx1 = (xl-xc);
	dx2 = (xr-xc);
	dy1 = (yl-yc);
	dy2 = (yr-yc);

	% TODO use dot here
	if (~issym(dx1))
		cosa = (dx1.*dx2+dy1.*dy2)./(hypot(dx1,dy1).*hypot(dx2,dy2));
	else
		h1 = sqrt(dx1.^2+dy1.^2);
		h2 = sqrt(dx2.^2+dy2.^2);
		cosa = (dx1.*dx2+dy1.*dy2)./(h1.*h2);
	end
end


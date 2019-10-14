% rough conversion to utm 
function [x,y,zone] = latlon2utm_simple(lat,lon)
	% point scale factor
	k0 = 0.9996;
	% number of zones
	nzone = 60;
	% degrees per zone
	dlon  = 360/nzone;
	% circumference at equator (scale by sin(x) for latitude
	% TODO
	%c    = 4.0075e7;
	c    = 6.378137e6*2*pi
	dx   = c/nzone;
	x0   = 5e5;
	l    = lon/dlon + 1/2*nzone;
	zone = floor(l)+1
	% the origin of the zones are centered, it at multiples of 3deg
	% and have an offset of 500km to avoid negative coordiantes
	x    = k0*(cosd(lat).*(frac(l)-1/2)*dx + x0);
	y0   = 0;
	% TODO zone letter
	%dlat  = 180/8;
	% northing is approximately the latitude scaled to m
	y    = k0*lat*c/360
	y    = ((y)) + y0
end


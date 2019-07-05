% Mon 26 Sep 13:34:31 CEST 2016
%
%% distance between corner of a triangle and its opposing mid-point
function [d2, d] = tri_distance_opposit_midpoint(X,Y)
	Xc = 0.5*(left(X)+right(X));
	Yc = 0.5*(left(Y)+right(Y));
	dX = X-Xc;
	dY = Y-Yc;
	d2 = dX.^2 + dY.^2;
	if (nargout() > 1)
		d = sqrt(d2);
	end
end


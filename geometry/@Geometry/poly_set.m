% Wed 10 Aug 15:58:49 CEST 2016
% Karl Kastner, Berlin
%% associate point at arbitary location with a polygon it is contained in
%% and assign the value of the polygon to it

% TODO, this should just return the index
% TODO arbitrary function (addition)
function [val, flag] = poly_set(poly,pval,x0,y0,val)
		if (nargin() < 5)
			val = NaN(size(x0));
		end
		flag = false(size(x0));

		% put polygon values to containing points,
		% last polygon has highest priority
		for idx=1:length(poly)
			X = poly(idx).X;
			Y = poly(idx).Y;
			if (X(1) ~= X(end) && Y(1) ~= Y(end))
				X(end+1) = X(1);
				Y(end+1) = Y(1);
			end
			fdx      = isfinite(poly(idx).X) & isfinite(poly(idx).Y);
			%fdx     = Geometry.contains(X(fdx),Y(fdx),x0,y0);
			fdx      = Geometry.inPolygon(X(fdx),Y(fdx),x0,y0);
			val(fdx) = pval(idx);
			flag(fdx) = true;
		end
end


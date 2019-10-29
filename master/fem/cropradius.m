function [x1 x2 y1 y2] = cropradius(P, v, tol)
	idx = find((v.^2 > tol^2));
	x1 = min(P(idx,1));
	x2 = max(P(idx,1));
	y1 = min(P(idx,2));
	y2 = max(P(idx,2));
end


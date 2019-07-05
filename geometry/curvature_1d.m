% 2018-08-17 19:01:58.506634041 +0200 geometry/curvature_1d.m
%% curvature of a sampled parametric curve in two dimensions
function c = curvature1d(x,y)
	if (~issym(x))
		dx     = cdiff(x);
		dy     = cdiff(y);
		d2y    = cdiffb(y,2);
		c = d2y ./ (dx.^2+dy.^2).^(3/2);
	else
		dy_dx  = diff(y,x);
		d2y_dx2 = diff(y,x,2);
		c = d2y_dx2 ./ (1+dy_dx.^2).^(3/2);
	end
end


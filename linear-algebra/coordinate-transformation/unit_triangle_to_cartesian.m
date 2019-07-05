% 2016-09-23 21:11:21.441722006 +0200
%% transform coordinates in unit triangle to cartesian coordinates
function xy_ = toUnitT(xy0,xy)
	% shift, rotate, scale, shear
	xyU = [0 1 0;
               0 0 1
               1 1 1];
	xy0 = [xy0; [1 1 1]];
	T = xy0 \ xyU;
	xy = [xy;ones(1,size(xy,2))];
	xy_ = xy*T;
	xy_=xy_(1:2,:);
end


% 2016-09-23 21:15:07.947414573 +0200
%% transform coodinates into unit triangle
function xy_ = toUnitT(xy0,xy)
	% shift, rotate, scale, shear
	xyU = [0 1 0;
               0 0 1
               1 1 1];
	xy0 = [xy0; [1 1 1]];
	T = xy0 \ xyU;
	xy = [xy;ones(1,size(xy,2))];
	xy_ = xy*inv(T);
	xy_=xy_(1:2,:);
end


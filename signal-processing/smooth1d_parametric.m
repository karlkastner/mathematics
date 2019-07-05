% Thu 22 Nov 16:48:08 CET 2018
%
%% smooth position of p0=x0,y0 between p1=x1,y1 and p2=x2,y2,
%% so that distance to p1 and p2 becomes equal
%% and the chord length remains the same
function [x0,y0] = smooth1d_parametric(x0,y0,x1,y1,x2,y2)
	load('mat/optimum-chord-length.mat','fun');
	l   = hypot(x1-x0,y1-y0) + hypot(x2-x0,y2-y0);
	x0_ = fun.x0(l,x1,x2,y1,y2);
	y0_ = fun.y0(l,x1,x2,y1,y2);

	% there are two solutions, choose the solution that is closest
	% to the original point, the other solution is mirrored
	d1  = hypot(x0-x0_(1,:),y0-y0_(1,:));
	d2  = hypot(x0-x0_(2,:),y0-y0_(2,:));
	fdx = d1<d2;
	x0(fdx)  = x0_(1,fdx);
	x0(~fdx) = x0_(2,~fdx);
	y0(fdx)  = y0_(1,fdx);
	y0(~fdx) = y0_(2,~fdx);
end


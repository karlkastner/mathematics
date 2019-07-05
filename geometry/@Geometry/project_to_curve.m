% Tue 27 Nov 10:37:00 CET 2018
%% closest point on a curve with respect to a point at distance to the curve
% TODO rename plumb line to project_to_segment
% TODO pass arguments to plumb line as joined vectors
function [xyp,id,p,d] = project_to_curve(xyc,xy0)
	% get closest mid point along curve
	id = knnsearch(mid(xyc),xy0);
	% for zero order, just take xyc(id,:)

	% project
	% make convex, i.e. limit to 0, 1, when closest point is on a protruding corner (or end point)
	cflag = true;
	[xp, yp, p, d]  = Geometry.plumb_line(xyc(id,1),xyc(id,2),xyc(id+1,1),xyc(id+1,2),xy0(:,1),xy0(:,2),cflag);
	xyp = [xp,yp];
end


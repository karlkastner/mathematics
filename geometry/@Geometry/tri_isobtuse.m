% Mo 2. Nov 18:11:55 CET 2015
% Karl Kastner, Berlin
%
%% flag obntuse triangles
% X : [x1,x2,x3] x-coordinates of triangle corners
% Y : [x1,x2,x3] y-coordinates of triangle corners
%
% TODO what with right triangles?
function [isobtuse Xc Yc] = tri_isobtuse(X,Y)
	[Xc Yc]  = Geometry.tri_excircle(X,Y);
	isobtuse = ~Geometry.inTriangle(X,Y,Xc,Yc);
end


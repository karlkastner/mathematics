% Di 3. Nov 15:35:22 CET 2015
% Karl Kastner, Berlin
%% smallest enclosing equilateral triangle with bottom site paralle to X axis
function [X Y elem s] = enclosing_triangle(X,Y)
	miX = min(X);
	maX = max(X);
	miY = min(Y);
	maY = max(Y);
	
	% centre of bottom side
	cX = 0.5*(miX+maX);
	cY = miY;

	% centre of left side	
	%[iX iY] = intersect([miX, maX], [maY maY],[cX,cX-cos(deg2rad(60))],[cY,cY+sin(deg2rad(60))]);
	[iX iY] = Geometry.intersect([miX, miX+cos(deg2rad(60))], [maY maY+sin(deg2rad(60))],[cX,cX-cos(deg2rad(60))],[cY,cY+sin(deg2rad(60))]);
%	d = [hypot(cX-miX,cY-maY);
%            maX-miX];
	% half side length
%	a = max(d);
	a = hypot(cX-iX,cY-iY)

	% the equilateral triangle
	X = [cX-a,cX,cX+a]';
	Y = [miY,miY+sqrt(3)*a,miY]';
	elem = [1 2 3]
	s = 2*a;
end


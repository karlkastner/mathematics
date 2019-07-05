% Di 3. Nov 15:43:57 CET 2015
% Karl Kastner, Berlin
%
%% intersect between two lines
function [px py] = intersect(x1,y1,x2,y2)
	b = [ x1(1)-x2(1); y1(1)-y2(1)];
	A = [(x1(2)-x1(1)),-(x2(2)-x2(1));
             (y1(2)-y1(1)),-(y2(2)-y2(1))];
	s = A \ b;
	px = x1(1)-s(1)*(x1(2)-x1(1));
	py = y1(1)-s(1)*(y1(2)-y1(1));
end


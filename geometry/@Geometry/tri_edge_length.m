% 2016-10-25 17:15:22.800900999 +0200
%
%% edge length of a triangle
function l = tri_edge_length(x,y)
	dx = right(x)-left(x);
	dy = right(y)-left(y);
	l = hypot(dx,dy);
end


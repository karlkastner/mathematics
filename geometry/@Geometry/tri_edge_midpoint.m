% Mon  5 Dec 12:27:24 CET 2016
%
%% mid point of a triangle
function [Xc Yc] = tri_edge_midpoint(X,Y)
	Xc = 0.5*(left(X)+right(X));
	Yc = 0.5*(left(Y)+right(Y));
end

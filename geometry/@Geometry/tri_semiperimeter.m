% 2016-10-25 17:23:53.619050012 +0200

%% semiperimeter of a triangle

function [s l] = tri_semiperimeter(x,y);
	l = Geometry.tri_edge_length(x,y);
	s = 0.5*sum(l,2);
end

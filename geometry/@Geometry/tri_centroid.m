% 2016-09-29 21:14:45.569318064 +0200
%
%% centroid of a triangle
function [xc yc] = tri_centroid(x,y)
	xc = mean(x,2);
	yc = mean(y,2);
end


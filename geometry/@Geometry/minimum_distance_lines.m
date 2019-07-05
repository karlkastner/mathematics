% Thu 27 Oct 13:31:31 CEST 2016
% Karl Kastner, Berlin
%
%% minimum distance of two lines in three dimensions
% point minimising sum of squared distances to a set of lines
%
%	parameter: xp,yp, s_i
%	min sum( (x_p - x_0i - s_i dx_i)^2 + (y_p - y_0i - s_i dy_i)^2 )
%	0 = df/dxp =   2 sum ( x_p - x_0i - s_i dx_i ) 
%	0 = df/dyp =   2 sum ( y_p - y_0i - s_i dy_i )
%	0 = df/dsi = - 2 dxi ( x_p - x_0i - s_i dxi) - 2 dyi (y_p - y_0i - s_i dyi)
%	           = - 2 dxi ( x_p - x_0i ) - 2 dyi (y_p - y_0i) + 1 s_i
%
function [xyp] = minimum_distance_lines(xy0,dir)
	n = size(dir,2);
	% factor 2 drops out
	A = [          n,          0, -dir(1,:);
                       0,          n, -dir(2,:);
	      -dir(1,:).', -dir(2,:).',    eye(n)];

	b = [ sum(xy0(1,:));
              sum(xy0(2,:));
	     -sum(dir.*xy0).' ];
	c = A \ b;
	xyp = c(1:2);
end


% Fri Feb 24 18:40:52 MSK 2012
% Karl Kästner, Berlin
% c.f. Dunavant
%
% baricentric coordinates and weights for gauss quadrature on the triangle
% midpoint rule
% 1 points, 2nd order accurate
function [w, b, flag] = int_2d_gauss_1()
		% integration weights
		w = 1;

		% integration points in barycentric coordiantes
		b = 1/3*[ 1 1 1 ];

		% mass matrix will be diagonal
		flag = 1;
end % int_2d_gauss_1


% Fri Feb 24 18:40:52 MSK 2012
% Karl Kästner, Berlin

% weigts for traiangle centre points scheme 2nd order accurate
% midpoint rule
function [w b flag] = int_2d_gauss_1()
		% integration weights
		w = 1;

		% integration points in barycentric coordiantes
		b = 1/3*[ 1 1 1 ];

		% mass matrix will not be diagonal
		flag = 0;
end % int_2d_tmp


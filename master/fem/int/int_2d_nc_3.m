% Fri Feb 24 18:41:11 MSK 2012
% Karl Kästner, Berlin

% trapezoidal rule for the triangle
function [w b flag] = int_2d_nc_3()
		% integration weights
		w = 1/3*[ 1 1 1 ]';

		% integration points in barycentric coordiantes
		b =  [ 1 0 0;
		       0 1 0;
		       0 0 1];

		% mass matrix will be diagonal
		flag = 1;
end % int_2d_nc_3()


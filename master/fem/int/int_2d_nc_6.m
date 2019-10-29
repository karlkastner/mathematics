% Fri Feb 24 18:41:11 MSK 2012
% Karl Kästner, Berlin

% weigts for triangle side mid point scheme, 2nd order accurate
function [w b flag] = int_2d_nc_6()
		% integration weights
		w = 1/3*[ 1 1 1 ]';

		% integration points in barycentric coordiantes
		b = 0.5*[ 0 1 1;
			  1 0 1;
			  1 1 0];

		% mass matrix will not be diagonal
		flag = 0;
end % int_2d_nc_6


% Fri Feb 24 18:41:11 MSK 2012
% Karl Kästner, Berlin

% weights for 3rd order accurate gauss quadrature
function [w b flag] = int_2d_gauss_3()
		% integration weights
		w = 1/3*[ 1 1 1]';

		% integration points in barycentric coordiantes
		b = 1/6*[ 4 1 1;
			  1 4 1;
			  1 1 4];
		% mass matrix will not be diagonal
		flag = 0;
end % int_2d_gauss_3


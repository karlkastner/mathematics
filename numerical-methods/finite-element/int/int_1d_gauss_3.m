% Wed Jul 11 17:27:33 MSK 2012
% Karl Kästner, Berlin

% coordinates and weights for numerical Gauss quadrature
% 5th-order accurate
function [w, b, flag] = int_1d_gauss_3()
	 %w = [0.888888888888889; 0.555555555555556];
	 %b = [0.000000000000000; 0.774596669241483];
         w = [ 0.444444444444444; 0.277777777777778; 0.277777777777778];
   	 b = [ 0.500000000000000   0.500000000000000;
   	       0.112701665379259   0.887298334620741;
               0.887298334620741   0.112701665379259];
	flag = 0;
end % int_1d_gauss_3


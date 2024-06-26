% Fri Jul 13 14:55:01 MSK 2012
% Karl Kästner, Berlin

% yields diagonal matrix for 35-point integration rule (24 zero coefficients), but not full accuracy
% 
function [w b flag] = int_3d_nc_11()
	w = 1/60*[1 1 1 1 4 4 4 4 4 4 32]';
	b = 1/4*[  4  0  0  0
                   0  4  0  0
                   0  0  4  0
                   0  0  0  4
                   0  0  2  2
                   0  2  0  2
                   2  0  0  2
                   0  2  2  0
                   2  0  2  0
                   2  2  0  0
                   1  1  1  1];
	flag = 11;
end


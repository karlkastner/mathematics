% Fri Jul 13 14:42:22 MSK 2012
% Karl Kästner, Berlin

% newton cotes rule
% yield diagonal mass matrix in combination with the 20-point tetrahedron (12 zero coefficients)
function [w b flag] = int_3d_nc_8()
	w = 1/40*[1; 1; 1; 1; 9; 9; 9; 9];
	b = 1/3*[3 0 0 0
                 0 3 0 0
                 0 0 3 0
                 0 0 0 3
                 0 1 1 1
                 1 0 1 1
                 1 1 0 1
                 1 1 1 0];
	flag = 20;
end


% Wed Jul 11 22:44:53 MSK 2012
% Karl Kästner, Berlin

% trapezoidal rule
% diagonal mass matrix in combination with the 4-point standard tetrahedron
function [w b flag] = int_3d_nc_4()
	w = [0.25; 0.25; 0.25; 0.25];
	b = [1 0 0 0;
             0 1 0 0;
             0 0 1 0;
	     0 0 0 1];
	% diagonal mass matrix scheme
	flag = 4;
end


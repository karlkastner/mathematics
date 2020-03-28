% Fri Feb 24 18:40:31 MSK 2012
% Karl Kästner, Berlin

% trapezoidal rule
function [w b flag] = int_1d_nc_2()
	% weights of integration points
	w = [0.5;
             0.5];
	% baricentric coordinates of integration points
	b = [1 0;
             0 1];
	% mark scheme to be diagonal
	flag = 1;
end % int_1d_nc_2()


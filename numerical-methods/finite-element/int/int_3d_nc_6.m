% Fri Feb 24 18:41:11 MSK 2012
% Karl Kästner, Berlin

% side mid-point integration rule
% newton cotes
% diagonal mass matrix in combination with the lagrangian 10-point tetrahedron
function [w b flag] = int_3d_nc_6()
	flag = 1;
	w = 1/6*ones(6,1);
	b = 0.5*[
	0     0     1     1
	0     1     0     1
	1     0     0     1
	0     1     1     0
	1     0     1     0
	1     1     0     0
	];
end % int_3d_nc_6


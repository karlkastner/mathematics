% Mon Jul 23 19:44:30 MSK 2012
% Karl Kästner, Berlin

% Simpson's rule (Kepplersche Faßregel)
function [w b flag] = int_1d_nc_3()
	b = [1.0 0.0;
             0.5 0.5
             0.0 1.0];
	w = 1/6*[1.0, 4.0, 1.0]';
	flag = 1;
end


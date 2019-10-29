% Mon Jul 23 19:48:41 MSK 2012
% Karl Kästner, Berlin


function [w b flag] = int_1d_nc_4()
	b = 1/3*[3 0;
                 2 1
                 1 2
		 0 3];
	w = 0.125*[ 1.0, 3.0, 3.0, 1.0]';
	flag = 1;
end % int_1d_nc_4
 


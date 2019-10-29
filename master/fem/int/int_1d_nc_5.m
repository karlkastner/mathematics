% Mon Jul 23 19:54:49 MSK 2012
% Karl Kästner, Berlin

function [w b flag] = int_1d_nc_5()
	b = 0.25*[4 0;
                  3 1;
                  2 2;
                  1 3;
                  0 4];
	w = 1/90*[ 7, 32, 12, 32, 7 ]';
	flag = 1;
end


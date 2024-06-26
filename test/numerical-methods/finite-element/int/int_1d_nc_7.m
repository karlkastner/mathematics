% Sun 11 Dec 00:40:11 CET 2016
% Karl Kästner, Berlin

function [w b flag] = int_1d_nc_7()
	b = 0.2*[6 0;
		 5 1;
                 4 2;
                 3 3;
                 2 4;
                 1 5;
                 0 6];
   	num =  [ 41 216  27  272 27 216 41 ]';
	den =  140;
	w   = num./den;

	flag = 1;
end % int_1d_nc_7()


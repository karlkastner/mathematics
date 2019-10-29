% Sun 11 Dec 00:40:11 CET 2016
% Karl Kästner, Berlin

% aka hardys rule
function [w b flag] = int_1d_nc_7()
	% coeffcients 2 and 5 are 0
	b = 0.2*[6 0;
		 5 1;
                 %4 2;
                 3 3;
                 %2 4;
                 1 5;
                 0 6];
%	w = 1/722*[95, 375, 250, 250, 375, 95]';
   	num =      [ 28 162 220 162 28 ]';
	den =  100*[ 1    1  1    1  1 ]';
	w   = num./den;

	flag = 1;
end % int_1d_nc_7()


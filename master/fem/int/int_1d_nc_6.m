% Mon Jul 23 19:57:49 MSK 2012
% Karl Kästner, Berlin

function [w b flag] = int_1d_nc_6()
	b = 0.2*[5 0;
                 4 1;
                 3 2;
                 2 3;
                 1 4;
                 0 5];
%	w = 1/722*[95, 375, 250, 250, 375, 95]';
   	num =   [ 19    25    25    25    25    19 ]';
	den =   [ 288    96   144   144    96   288 ]';
	w   = num./den;

	flag = 1;
end % int_1d_nc_6()


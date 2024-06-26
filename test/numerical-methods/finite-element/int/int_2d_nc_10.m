% Mon Jul 23 19:32:35 MSK 2012
% Karl Kästner, Berlin

function [w b flag] = int_2d_nc_10()
	b = 1/3*[
	    3    0    0
	    2    1    0
	    1    2    0
	    0    3    0
	    2    0    1
	    1    1    1
	    0    2    1
	    1    0    2
	    0    1    2
	    0    0    3 ];
	num = 2*[ 1     3     3     1     3     9     3     3     3     1 ]';
	den = [ 60    80    80    60    80    40    80    80    80    60 ]';
	w = num./den;
	flag = 1;
end % int_2d_nc_10()


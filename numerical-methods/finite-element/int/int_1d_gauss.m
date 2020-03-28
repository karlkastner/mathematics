% Do 18. Jun 16:34:49 CEST 2015
% Karl Kastner, Berlin
% function [w, b, flag] = int_gauss_1d(order)
function [w, b, flag] = int_gauss_1d(order)
	switch(order)
		case {1}
			[w, b, flag] = int_1d_gauss_1();
		case {2}
			[w, b, flag] = int_1d_gauss_2();
		case {3}
			[w, b, flag] = int_1d_gauss_3();
		case {4}
			[w, b, flag] = int_1d_gauss_4();
		case {5}
			[w, b, flag] = int_1d_gauss_5();
		case {6}
			[w, b, flag] = int_1d_gauss_6();
		otherwise
			error('int_gauss_1d');
	end
end % int_gauss_1d


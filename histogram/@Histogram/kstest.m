% 2014-11-22 21:35:49.869311850 +0100
% Karl Kastner, Berlin

function [alpha obj] = kstest(obj,h1,h2,n1,n2)
	H1 = cumsum(h1);
	H2 = cumsum(h2);
	d = max(abs(H1-H2));
	k = d/sqrt((n1+n2)/(n1*n2));
	alpha = 1-kolmcdf(k);
end


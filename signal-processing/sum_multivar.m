% 2015-08-14 14:57:06.895344109 +0200
%% sum of matrix entries of bivariate ar1 process
function s = sum_multivar(r1,r2,m)
	s = m*(1 + r1/(1-r1) + r2/(1-r2)) ...
		- r1*(1-r1^m)/(1-r1)^2 ...
		- r2*(1-r2^m)/(1-r2)^2;
end


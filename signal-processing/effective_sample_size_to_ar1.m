% Fri  2 Feb 13:45:40 CET 2018
%% convert effective sample size to ar1 correlation
function rho = effective_sample_size_to_ar1(n,twosided)
	if (~twosided)
		rho = (1-n)/(1+n);
	else
	rho = (n - 3)/(3*(n + 1)) ...
		+ (2^(1/3)*n^(1/3)*(3*3^(1/2)*(n^2 + 27)^(1/2) + 5*n^2 + 3*3^(1/2)*n*(n^2 + 27)^(1/2) + 27)^(1/3))/(3*(n + 1)) ...
		- (2^(2/3)*n^(2/3)*(n + 9))/(3*(n + 1)*(3*3^(1/2)*(n^2 + 27)^(1/2) + 5*n^2 + 3*3^(1/2)*n*(n^2 + 27)^(1/2) + 27)^(1/3));
	end
end


% 2015-08-07 10:25:23.125806143 +0200
%% variance of an autocorrelated finite process
function v = ar1_var_factor_1(rho,n,m)
%	N = 
	v = 1/m^2*sum_ii(rho,m) - 2/(n*m)*sum_ij(rho,n,m) + 1/n^2*sum_ii(rho,n);
end


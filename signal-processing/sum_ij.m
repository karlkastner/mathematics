% Fr 7. Aug 09:48:26 CEST 2015
% Karl Kastner, berlin
%% sum of ar1 matrix
%% sum_{i=1}^n sum_{j=1}^m r^|i-j|
function s = sum_ij(rho,n,m)
%	s = 1/(1-rho)*(m*(1+rho) - rho*(1-rho^m)/(1-rho) + - rho^(n+1)/(1-rho)*(rho^-m - 1));
	s = m*(1+rho)/(1-rho) - rho/(1-rho)^2*( 1 - rho^m + rho^(n-m) - rho^n);
end



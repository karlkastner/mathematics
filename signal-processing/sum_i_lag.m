% 2015-08-05 15:49:54.653409204 +0200
% Karl Kastner, Berlin

%% sum of ar1 matrix with lag
%% sum_i=1^n rho^|i-k|

function s = sum_i_lag(rho,n,k)
	s(1) = (1+rho)/(1-rho) - 1/(1-rho)*(rho^k+rho^(n-k+1));
end


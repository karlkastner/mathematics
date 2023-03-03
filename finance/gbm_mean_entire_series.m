% Fri 24 Feb 11:41:44 CET 2023
%
% average of the geometric brownian for the entire series, not just the final state
% E(S,0,T) = S0/T int_0^t exp(mu t) dt
%	 = S0/(mu T) exp(mu t) |_0^T
%	 = S0/(mu T) (exp(mu T) - 1) 
function mu_S = gbm_mean_entire_series(S0,r,sigma,T)
	mu = r + 0.5*sigma^2;
	mu_S = S0./(mu*T)*(exp(mu*T) - 1);
end


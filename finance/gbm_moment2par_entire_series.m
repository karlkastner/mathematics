function [r,sigma] = gbm_moment2par_entire_series(S0,mu_S,sd_S,T)

%	mu_S = S0./(mu*T)*(exp(mu*T) - 1);
%	1+mu_S*mu*T/S0  = exp(mu*T)
%	(1+mu_S*mu*T/S0)^2  = exp(2*mu*T)

%	s2_S = S0^2/T*( 1/(2*mu+sigma^2)*(exp(2*mu*T + sigma^2*T) - 1) - 1/(2*mu)*(exp(2*mu*T)-1));
%	s2_S = S0^2/T*( 1/(2*mu+sigma^2)*(exp(2*mu*T)*exp(sigma^2*T) - 1) - 1/(2*mu)*(exp(2*mu*T)-1));
%syms S0 T mu sigma mu_S s2_S
%	eq = simplify(s2_S - S0^2/T*( 1/(2*mu+sigma^2)*((1+mu_S*mu*T/S0)^2*exp(sigma^2*T) - 1) - 1/(2*mu)*((1+mu_S*mu*T/S0)^2-1)))
	resfun  = @(par) [(mu_S - gbm_mean_entire_series(S0,par(1),par(2),T)); (sd_S - gbm_std_entire_series(S0,par(1),par(2),T))];
	[r0,s0] = gbm_moment2par(S0,mu_S,sd_S,0.5*T);
	par = lsqnonlin(resfun,[r0,s0]);
	r = par(1);
	sigma = par(2);
end


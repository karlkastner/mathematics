% Mon 23 Jan 14:17:04 CET 2023
%
% transform moments of the geometric brownian motion to its parameters
% input:
%	S0    : S(0)
%	t     : time at the moments of S(t) were determined
%	mu_St : E[S(t)]	mean of S(t)
%	sd_St : std(S(t))	standard deviation of S(t)
% output:
%	r     : risk free interest rate
%	sigma : volatility
function [r,sigma] = gbm_moment2par(mu_St,sd_St,S0,t)
	% mu_s = S0*exp(mu*t)
	% sd_s = S0*exp(mu*t).*sqrt(exp(sigma.^2*t)  - 1);
	% sd_s/mu_s = sqrt(exp(s^2*t) - 1))
	% log(1+(sd_s/mu_s)^2) = s^2*t
	% s^2*t = (log(1+sd^2/mu_s^2));
	sigma2 = log(1+(sd_St/mu_St)^2)./t;
	sigma  = sqrt(sigma2);
	mu     = log(mu_St/S0)./t;
	r      = mu - 0.5*sigma.^2;
end


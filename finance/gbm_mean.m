% Sun 12 Jan 10:47:33 +08 2020
% mean of the geometric brownian motion with drift, not equal to median
% function s_mu = gbm_mean(t,r,sigma,S0)
% r : risk free interest rate
function s_mu = gbm_mean(t,r,sigma,S0)
	% perceived interest rate due to the presence of risk
	mu = r + 1/2*sigma^2;
	s_mu = S0*exp(mu*t);
end


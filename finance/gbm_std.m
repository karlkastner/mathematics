% Sun 12 Jan 10:47:33 +08 2020
% standard deviation of the geometric brownian motion with drift, not equal to median
function s_std = gbm_mean(t,r,sigma,S0)
	mu = r+1/2*sigma.^2;
	s_std = S0*exp(mu*t).*sqrt(exp(sigma.^2*t)  - 1);
end


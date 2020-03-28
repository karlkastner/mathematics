% Sun 12 Jan 11:22:39 +08 2020
% function f = gbm_cdf(t,s,r,sigma,S0)
% r : risk free interest rate
function f = gbm_cdf(t,s,r,sigma,S0)
	%mu_    = log(S0) + (r-1/2*sigma.^2).*t;
	mu     = r + 0*1/2*sigma.^2;
	mu_    = log(S0) + mu.*t;
	sigma_ = sigma*sqrt(t);
	f      = logncdf(s,mu_,sigma_);
end


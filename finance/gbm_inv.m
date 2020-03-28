% Sun 12 Jan 11:23:55 +08 2020
% function f = gbm_inv(t,p,r,sigma,S0)
% r : risk free interest rate
function f = gbm_inv(t,p,r,sigma,S0)
	mu      = r       + 1/2*sigma.^2;
	S_mu    = log(S0) + (r-1/2*sigma.^2).*t;
	S_sigma = sigma*sqrt(t);
	f       = logninv(p,S_mu,S_sigma);
end


% Sun 12 Jan 10:52:26 +08 2020
% median of the geometric brownian motion with drift, not equal to median
% r : risk free interest rate
function s_mu = gbm_median(t,r,sigma,S0)
	% r = mu-1/2*sigma^2
	s_mu = S0*exp(r*t);
end


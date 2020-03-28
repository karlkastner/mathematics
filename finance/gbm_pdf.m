% Sun 12 Jan 10:42:17 +08 2020
% Karl Kastner, Berlin
%
% probability density of geometric browniam motion with drift
%
% t     : time for which pdf to compute
% S     : stock price for which probability to compute
% r     : risk free interest rate, the amount at which the price increases or falls
%	  when there is no risk (sigma=0)
% sigma : standard deviation in log-space
% S0    : stock price at t = 0
function f = gbm_pdf(t,S,r,sigma,S0)
%	f  = 1./(sqrt(2*pi*t)*s*sigma)*exp( -(log(s) - (log(S0)+(r-1/2*sigma^2)*t)).^2./(2*sigma^2*t) );
% lognormal :
	%mu_    = log(S0) + (r-1/2*sigma.^2).*t;
	mu     = r+0*1/2*sigma.^2;
	mu_    = log(S0) + mu.*t;
	sigma_ = sigma.*sqrt(t);
%	x      = s
%	f(:,2)  = 1./(sqrt(2*pi)*x.*sigma_).*exp( -(log(x) - mu_).^2./(2*sigma_^2) );
	f      = lognpdf(S,mu_,sigma_);
end


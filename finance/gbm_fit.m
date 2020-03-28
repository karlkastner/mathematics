% Sun 12 Jan 11:46:52 +08 2020
%
% [r,sigma,S0,serr_] = gbm_fit(t,S,mode)
%
% fit parameters of the geometric brownian motion (black-scholes) model
%
% input:
% t : vectors of sample times
% S : vector of sample prices
%
% output:
% r  : interest rate, not risk free
%      perceived interest rate is mu = r+1/2*sigma^2
% sgima1  : volatility for unit time step
% S0 : fitted price at t0
%
% sigma_t^2 = S_0^2 exp(2 mu t) (exp(sigma^2t) - 1)
%           = S_t^2 (exp(-sigma^2 t) - 1)
% log(1 + res^2/S_t^2) = sigma^2 t
%
%function [r,sigma1,S0,serr_] = gbm_fit(t,S)
function [fit] = gbm_fit(t,S)
	nt = length(t);
	lS = log(S);

		% c.f. Teka, 2012
		% difference of increments of logarithms is are normally and independently distributed
		% dS = r*S*dt + random part
		%r  = S(1:end-1) \ diff(S); % \ S(1:end-1);
		% S_i+1 = S_i exp(r*dt)
		%r  = mean(log(S(2:end)./S(1:end-1))); % S0 = S(1);
		%r  = log(S(2:end)./S(1:end-1))); % S0 = S(1);
		dt    = diff(t);
		n     = length(t);
		lrat     = log(S(2:end)./S(1:end-1))./dt;
		mu       = mean(lrat);
		s2       = mean(((lrat-mu).^2)./dt);
		sigma1   = sqrt(s2);
		skew     = skewness(lrat);
		kurt     = kurtosis(lrat);
		% the risk free interest rate
		% this correction does no seem necessary here, as the single step
		% returns have zero mean influence by volatility
		%r     = mu-1/2*sigma^2;
		r     = mu;
		serr_ = [sigma1/sqrt(n-1), sigma1/sqrt(2*n)];
		% sigma is xi^2 with distribution n/xi^2(n,1-p) sigma < sigma < n/xi^2(n,1-p) sigma
		% the chi^2 distribution has variance 2 n, so variance of sigma is sigma/sqrt(2n) sigma
		%S0 = S(1);
		S0    = mean(S.*exp(mu*(t(1)-t)));

	Sp = S0*exp(r*(t-t(1)));
	res = Sp-S;

	fit = struct();
	fit.rate = r;
	fit.sigma = sigma1;
	fit.skewness = skew;
	fit.kurtosis = kurt;
	fit.serr_ = serr_;
	fit.S0 = S0;
end % gbm_fit

	

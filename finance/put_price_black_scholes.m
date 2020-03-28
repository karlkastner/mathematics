% Sun 19 Jan 14:44:42 +08 2020
% P : price for option to put (sell) stock at time T for price K
% S : stock price at the moment the option is purchased
% r : interest rate (risk free)
% s : volatility
% function P = put_price_black_scholes(K,S0,r,s,T)
function [P,C] = put_price_black_scholes(K,S0,r,s,T)
	d1 = 1./(s.*sqrt(T)).*(log(S0./K) + (r+1/2*s.^2).*T);
	d2 = d1 - s*sqrt(T);

%d2_ = d2;
%d2 = d1;
%d1 = d2_;

% probability of exercising vs probability of not exercising

% Note : both d1 and d2 are probabilities that the option is exectuded (same sign)
% not as often claimed d1 for execution and d2 for not execution,
% the only difference is the risk neutrality (d2 is risk neutral)
%	but why is it then r-1/2 sigma^2, not simply r?



%r_ = r+0*1/2*s^2;
exp(r*T)*S0*normcdf(d1)
%exp(r_*T)*S0*normcdf(d2)
% E[S|S>K]*P(S>K)


%quad(@(S_T) S_T.*gbm_pdf(T,S_T,r,s,S0),K,K*10)./(1-gbm_cdf(T,K,r,s,S0))
%quad(@(S_T) S_T.*gbm_pdf(T,S_T,r,s,S0),K,K*10) %./(1-gbm_cdf(T,K,r,s,S0)).*(1-gbm_cdf(T,K,r,s,S0))
quad(@(S_T) S_T.*gbm_pdf(T,S_T,r-1/2*s^2,s,S0),K,K*10) %./(1-gbm_cdf(T,K,r,s,S0)).*(1-gbm_cdf(T,K,r,s,S0))
%quad(@(S_T) S_T.*gbm_pdf(T,S_T,r,s,S0),sqrt(eps)*K,K) %./(gbm_cdf(T,K,r,s,S0)).*(gbm_cdf(T,K,r,s,S0))

normcdf(d2)
quad(@(S_T) gbm_pdf(T,S_T,r-1/2*s^2,s,S0),K,K*10) %./(1-gbm_cdf(T,K,r,s,S0)).*(1-gbm_cdf(T,K,r,s,S0))

pause

% Nd2 : risk neutral probability, that option is exectuted


%	gbm_cdf(T,K,r,s,S0)
%	normcdf(-d1)
%	quad(@(K) K.*gbm_pdf(T,K,r,s,S0),1e-4*K,K)
%./(gbm_cdf(T,K,r,s,S0))
	
%	1-gbm_cdf(T,K,r-1/2*s^2,s,S0)
%	normcdf(-d1)
%	1-normcdf(-d1)
%	normcdf(-d2)
%pause
%K

%	normcdf(-d2).*K
%	normcdf(d2).*K
%	quad(@(K) K.*gbm_pdf(T,K,r,s,S0),K,10)
	%	*/(1-gbm_cdf(T,K,r,s,S0))
%pause
%
	% strike price is discounted by r, as it incurrs no risk
	%				  probability, that S_T<K, i.e. option is executed
	%				  P(S_T<K)*S0
	P = normcdf(-d2).*K.*exp(-r.*T) - normcdf(-d1).*S0;
	C = normcdf(+d1).*S0            - normcdf(+d2).*K.*exp(-r*T);
% call+put : P+C = -K*exp(-r+T) - S0, with N(x) = 1-N(-x)
end


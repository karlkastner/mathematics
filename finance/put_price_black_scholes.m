% Sun 19 Jan 14:44:42 +08 2020
% P : price for option to put  (sell) stock at time T for price K
% C : price for option to call (buy) stock at time T for price K
% S : stock price at the moment the option is purchased
% r : interest rate (risk free)
% s : volatility
% function P = put_price_black_scholes(K,S0,r,s,T)
function [P,C] = put_price_black_scholes(K,S0,r,s,T)
	d1 = 1./(s.*sqrt(T)).*(log(S0./K) + (r+0.5*s.^2).*T);
	d2 = d1 - s*sqrt(T);

	% call+put : P+C = -K*exp(-r+T) - S0, with N(x) = 1-N(-x)
	P = normcdf(-d2).*K.*exp(-r.*T) - normcdf(-d1).*S0;
	C = normcdf(+d1).*S0            - normcdf(+d2).*K.*exp(-r*T);
end


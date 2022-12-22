% Mon  4 Jul 10:58:48 CEST 2022
% x : start value
% b : drift rate
% a : threshold to pass at t=tau
%
% c.f. Abundo
function P = brownian_drift_hitting_probability(t,a,y0,b,s)
% ends below threshold
% ends above threshold but was below
	% P(tau) <= t
	P = 1 - normcdf((a-y0-b*t)./(s*sqrt(t))) - exp(-2*b*(a-y0)*s)*normcdf(-(a-y0-b*t)./(s*sqrt(t)));
%	P = 1 + normcdf(a,b*t+y0,s*sqrt(x)) - exp(2*b*(a-y0)*s)*normcdf((b*t - (a-y0))./(s*sqrt(t)));
end


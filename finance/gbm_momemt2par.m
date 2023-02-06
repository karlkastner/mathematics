% Mon 23 Jan 14:17:04 CET 2023
function [mu,sd] = gbm_moment2par(E,sqrt_V,t)
	% E = exp(mu*t)
	% V = exp(2*mu*t)*(exp(s^2*t) - 1)
	% V/E^2 = exp((s^2*t) - 1)
	% s^2*t = log(1-sd^2/E^2);
	s2 = log(1-(sqrt_V/E).^2)./t;
	sd  = sqrt(s2);
	mu = log(E)./t;
end


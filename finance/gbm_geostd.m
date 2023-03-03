% Sun 12 Jan 10:47:33 +08 2020
% standard deviation of the geometric brownian motion with drift, not equal to median
function s_std = gbm_geostd(t,r,sigma,S0)
	s_std = exp(sqrt(t)*sigma);
end


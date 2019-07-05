% Mi 2. Mär 14:48:22 CET 2016
%
%% transform parameters to mode (mu, sd) for the log normal distribution
function [mu, sd] = logn_param2mode(lmu,lsd)
	mu = exp(lmu + 0.5*lsd^2);
	sd = sqrt(exp(lsd^2)-1)*exp(lmu+0.5*lsd^2);
end


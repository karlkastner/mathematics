% 2016-10-23 20:54:10.409476622 +0200
% Karl Kastner, Berlin
%
%% transform modes (mu,sd) to parameters of the gamma distribution
%
function [a, b] = gamma_mode_to_parameter(mu,sd)
	% mean     = a*b
	% variance = sd^2 = a*b^2
	b = sd^2/mu;
	a = mu/b;
end


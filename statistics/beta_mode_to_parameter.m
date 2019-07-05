% 2016-10-23 21:12:46.030078012 +0200
% Karl Kastner, Berlin
%
%% transform modes (mean and sd) to paramets of the beta function
function [a b] = beta_mode_to_parameter(mu,sd)
	s2 = sd*sd;
	a = mu*(mu*(1-mu)/s2-1);
	b = (1 + (mu*mu - mu)/s2)*(mu - 1);
end


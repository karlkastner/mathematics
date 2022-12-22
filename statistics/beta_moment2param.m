% 2016-10-23 21:12:46.030078012 +0200
% Karl Kastner, Berlin
%
%% transform central moments (mean and sd) to parameters of the beta function
function [a, b] = beta_moment_to_parameter(mu,sd)
	s2 = sd*sd;
	a = mu*(mu*(1-mu)/s2-1)
	b = (1 + (mu*mu - mu)/s2)*(mu - 1)

%	b = (mu^3 - 2*mu^2 + mu*sd^2 + mu - sd^2)/sd^2
%	a =  b*mu/(1 - mu)
end


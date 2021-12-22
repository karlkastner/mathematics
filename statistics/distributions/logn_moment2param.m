% Mi 2. Mär 14:47:56 CET 2016
% Karl Kastner, Berlin
%
%% transform modes (mu,sd) to parameters of the log normal distribution
function [lmu,lsd] = logn_moment2param(mu,sd)
	lsd = sqrt(log(1 + sd.^2./mu.^2));
	lmu = log(mu)-0.5*lsd.^2;
end

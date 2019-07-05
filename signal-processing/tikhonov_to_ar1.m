% Fri  2 Feb 13:17:48 CET 2018
%% convert coefficient of the tikhonov regularization to correlatioon of the ar1 process
function [rho, diag_]= titkhonov_to_ar1(lambda)
	%a = (2*lambda - (4*lambda + 1)^(1/2) + 1)/(2*lambda)
	rho   = 1 + (1 - (4*lambda + 1)^(1/2))/(2*lambda);
	diag_ = 1/(1-rho)^2*[-rho,1+rho^2,-rho];
end


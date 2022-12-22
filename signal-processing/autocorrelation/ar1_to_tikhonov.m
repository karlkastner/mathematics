% Fri  2 Feb 13:17:38 CET 2018
%% convert ar1 correlation to tikhonovs lambda
function [lambda] = ar1_to_tikhonov(rho)
	lambda = rho/(1-rho)^2;
end


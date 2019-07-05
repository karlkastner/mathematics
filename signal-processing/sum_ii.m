% 2015-06-22 17:39:23.854197350 +0200
% Karl Kastner, Berlin
%
%% sum of ar1 matrix
%% sum_i=1^n sum_j=1^n rho^|i-j|
%% this is for the variance, take square root for the standard deviation factor
%
%
function [s] = sum_ii(r,n)
	s = 1/(1-r)^2*(n*(1-r*r)  + 2*r*(r.^n - 1));
end


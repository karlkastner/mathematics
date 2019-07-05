% Sa 1. Aug 12:58:05 CEST 2015
% Karl Kastner, Berlin
%
%% correlate to correlated standard normally distributed vectors
function [x, y] = randc(rho,n)
	x = randn(n,1);
	y = rho*x + sqrt(1-rho^2)*randn(n,1);
end % randc

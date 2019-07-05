% Mo 15. Feb 11:49:20 CET 2016
% Karl Kastner, Berlin
%% autocorrelation of the ar2 process
%% X_i + a1 X_i-1 + a2 X_i-2 = 0
function [rho1,rho2] = ar2acf(a)
	rho(1) = a(1)/(1-a(2));
	rho(2) = (a(1)^2 + a(2) - a(2)^2)/(1-a(2));
	

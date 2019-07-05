% 2015-08-05 11:50:08 +0200
% Karl Kastner, Berlin
%
%% error variance of a single sample of a finite length ar1 process
%% with respect to the mean, averaged over the population
%
% function s2 = varar1(s,rho,n)
%
% Note for testing: use var(x,1), as this is the full population
function s2 = varar1(s,rho,n)
%	s2 = 1 - 1/n^2*1/(1-rho)^2*(1+rho-rho^(n-1)*(1+rho^2));
	s2 = s^2-mu2ar1(s,rho,n);
end


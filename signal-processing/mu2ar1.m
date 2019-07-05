% Mi 5. Aug 11:43:37 CEST 2015
% Karl Kastner, Berlin
%
%% error variance of the mean of the finite length ar1 process
%%
%%
%% (mu)^2 = (sum epsi)^2 = sum_i sum_j eps_i eps_j = sum_ii(rho,n)/n^2
%% this has the limit s^2 for rho->1 
%
%function mu2 = mu2ar1(s,rho,n)
function mu2 = mu2ar1(s,rho,n)
	mu2 = 1/n*s^2*((1+rho)/(1-rho) ...
                  -2/n*1/(1-rho)^2*rho*(1-rho^n));
end


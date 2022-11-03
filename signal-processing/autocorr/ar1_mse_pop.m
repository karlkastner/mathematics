% Mi 7. Okt 14:45:45 CEST 2015
% Karl Kastner, Berlin
%
%% variance of the population mean of a single realisation around zero
%%
%% E[(mu_N-0)^2] = E[mu_N^]
%%
function s2 = ar1_var_pop(sigma,rho,N)
	s2 = sigma^2*(  1/N*(1+rho)/(1-rho) ...
		      - 2*rho/N^2*(1-rho^N)/(1-rho)^2);		
end


% Mi 7. Okt 15:02:47 CEST 2015
% Karl Kastner, Berlin
% 
%% standard error of a single sample of an ar1 correlated process
function s2 = ar1_var_mu_single_sample(sigma,rho,n,idx)
	s2 = sigma^2 ...
             - 2/n*mid_term_single_sample(sigma,rho,idx,n) ...
             + ar1_mse_pop(sigma,rho,n);
end



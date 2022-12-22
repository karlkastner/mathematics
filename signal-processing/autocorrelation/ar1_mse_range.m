% Di 13. Okt 14:51:57 CEST 2015
% Karl Kastner, Berlin
%
%% mean standard error of the mean of a range of values taken from an ar1 process
function s = ar1_mse_range(sigma,rho,n,m)
	s = sigma^2*ar1_var_factor(rho,n,m);
end


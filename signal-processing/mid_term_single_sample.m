% 2015-10-07 15:13:10.469295036 +0200
% Karl Kastner, Berlin
%% variance of single sample, mid term
function s2 = mid_term_single_sample(sigma,rho,idx,n)
	s2 = 	sigma^2 * ( ...
		1/(1-rho)*(rho*(1-rho.^(idx-1)) + (1-rho.^(n-idx+1))) );
end


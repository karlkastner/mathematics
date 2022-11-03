% Mi 7. Okt 14:40:59 CEST 2015
% Karl Kastner, Berlin

function S = test_ar1_var_mu_single_sample()

sigma = pi;
rho = 1/3;
N = 20;

S = [];
for idx=1:N
	S(idx,1) = mid_term_single_sample(sigma,rho,idx,N);	
	S(idx,2) = mid_term_single_sample_(sigma,rho,idx,N);
end
S

S = [];
for idx=1:N
	S(idx,1) = ar1_var_mu_single_sample(sigma,rho,idx,N);
	S(idx,2) = ar1_var_mu_single_sample_(sigma,rho,idx,N);
end
S
end % function

function s2 = ar1_var_mu_single_sample_(sigma,rho,idx,n)
	s2 = sigma^2 ...
             - 2/n*mid_term_single_sample_(sigma,rho,idx,n) ...
             + ar1_var_pop(sigma,rho,n); % no underscore here, tested in recursive test
end

% E[(ei - mu)^2] = E[(ei - 1/n sum e_j)^2]
%                = E[ei^2] - 2/n E[ei sum e_j] + E[(1/n sum e_i)^2]
%                = s^2  - 2/n mid_term + sp^2
function s = mid_term_single_sample_(sigma,rho,idx,n)
	N = 1:n;
	D = N-idx;
	R = rho.^abs(D);
	s = sigma^2*sum(R);
end


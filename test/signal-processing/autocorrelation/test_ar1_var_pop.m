% Mi 7. Okt 14:43:17 CEST 2015
% Karl Kastner, Berlin

function failed = test_ar1_var_pop()

N = 20;
rho = 1/3;
sigma = 1;

s2(1) = ar1_var_pop(sigma,rho,N)
s2(2) = ar1_var_pop_(sigma,rho,N)

failed = (s2(1) ~= s2(1))

end

% s_p^2 = E[mu^2] - E[mu]^2
%       = (1/N sum ei)(1/N sum ej)
function s2 = ar1_var_pop_(sigma,rho,N)
	I  = repmat((1:N)',1,N);
	J  = I';
	D  = I-J;
	R  = rho.^abs(D);
	s2 = sigma^2*sum(sum(R))/N^2;
end % ar1_var_pop


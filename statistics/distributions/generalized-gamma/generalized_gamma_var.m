function s2 = generalized_gamma_var(a,d,p)
	s2 = a^2*(gamma((d+2)/p)/gamma(d/p) - (gamma((d+1)/p)/gamma(d/p))^2);
end


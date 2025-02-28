function sd = generalized_gamma_std(a,d,p)
	s2 = generalized_gamma_var(a,d,p);
	sd = sqrt(sd);
end


function sd = powerpdf_std(a)
	s2 = powerpdf_var(a);
	sd = sqrt(s2);
end


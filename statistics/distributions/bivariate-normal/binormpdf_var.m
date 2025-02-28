function s2 = binormpdf_var(p,mu1,mu2,sd1,sd2)
	s2 = p*sd1^2 + (1-p)*sd2^2 + p*(1-p)*(mu1 - mu2)^2;
end


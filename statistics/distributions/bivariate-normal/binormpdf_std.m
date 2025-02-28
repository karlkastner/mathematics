% Tue  4 Feb 16:33:31 CET 2025
function sd = binormpdf_std(p,mu1,mu2,sd1,sd2)
	% \sigma^2 = p sigma1^2 + (1-p) sigma2^2 + p*(1-p)*(mu0 - mu1)^2
	s2 = binormpdf_var(p,mu1,mu2,sd1,sd2);
	sd = sqrt(s2);
end


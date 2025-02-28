% Tue  4 Feb 16:33:40 CET 2025
function mu = binormpdf_mean(p,mu1,mu2,sd1,sd2)
	mu = p*mu1 + (1-p).*mu2;
end

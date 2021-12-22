function sd = logn_std(lmu,lsd)
	mu = logn_mu(lmu,lsd);
	sd = sqrt(exp(lsd.^2)-1).*exp(lmu+0.5*lsd.^2);
end

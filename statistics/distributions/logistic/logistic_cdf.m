function c = logistic_cdf(x,mu,a)
	c = 1-1./(1 + exp(a*(x-mu)));
end


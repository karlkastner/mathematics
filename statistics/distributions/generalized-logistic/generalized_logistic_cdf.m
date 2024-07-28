function c = generalized_logistic_cdf(x,mu,a,b)
	c = 1-1./(1 + exp(a.*(x-mu))).^b;
end


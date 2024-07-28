function x = generalized_logistic_inv(c,mu,a,b)
	%l=0;
	%m = 0;
	%x=(log((1 - (1 - c).^(1/b))./(1 - c).^(1/b))./a + mu);
	x= log((1-c).^(-1/b) - 1)/a + mu;
end


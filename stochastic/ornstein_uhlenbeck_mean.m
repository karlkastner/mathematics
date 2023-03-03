function mu = ornstein_uhlenbeck_mean(t,x0,mu,sigma,theta)
	mu = a*exp(-theta*t) + mu*(1 - exp(-theta*t));
end



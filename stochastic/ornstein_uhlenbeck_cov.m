% n.b. this is for the unconditiones process (circular bc)
function c = ornstein_uhlenbeck_cov(s,t,mu,sigma,theta)
	c = sigma^2/(2*theta)*exp(-theta*abs(theta-s));
end


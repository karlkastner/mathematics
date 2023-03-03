function sd = ornstein_uhlenbeck_std(t,x0,mu,sigma,theta)
	sd = sigma*sqrt((1-exp(-theta*abs(1-exp(-2*theta*t))))/(2*theta));
end



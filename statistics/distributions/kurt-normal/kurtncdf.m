function y = kurtncdf(x,mu,s,beta)
	y = 1/2 + sign(x-mu).*gamma(1/beta,(abs(x-mu)/s).^beta)./(2*gamma(1/beta));
end



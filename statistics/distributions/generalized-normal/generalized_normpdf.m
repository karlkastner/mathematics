function p = generalized_normpdf(x,mu,a,b)
	p = b./(2*a.*gamma(1./b)).*exp(-(abs(x-mu)./a).^b);
end


function p = normpdf_sym(x,mu,s)
	p = 1./(sqrt(2*pi)*s)*exp(-(x-mu).^2/(2*s.^2));
end

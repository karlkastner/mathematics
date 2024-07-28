function c = generalized_normcdf(x,mu,a,b)
	%c = 1/2 + sign(x-mu).*1./(2*gamma(1./b).^0).*gammainc(abs((x-mu)./a).^b,1./b);
	c = 0.5*(1 + sign(x-mu).*gammainc(abs((x-mu)./a).^b,1./b));
end


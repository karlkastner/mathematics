function p = lognpdf_sym(x,a,b)
	if (issym(x) || issym(a) || issym(b))
		pi_ = sym('pi')
	else
		pi_ = pi;
	end
	p = 1./(x*b*sqrt(2*pi_)).*exp(-(log(x)-a).^2./(2*b.^2));
end

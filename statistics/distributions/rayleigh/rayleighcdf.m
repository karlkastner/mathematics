function p = rayleighcdf(x,s)
	p = x./(s.^2)*exp(-x.^2./(2*s.^2));
end


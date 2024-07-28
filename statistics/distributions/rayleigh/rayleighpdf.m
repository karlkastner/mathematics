% note that this is just a weibull dist
function p = rayleighpdf(x,s)
	p = x./(s.^2).*exp(-x.^2./(2*s.^2));
end


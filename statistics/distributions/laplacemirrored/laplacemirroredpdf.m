% Sat 29 Jun 10:29:30 CEST 2024
function p = laplacewrappedpdf(x,mu,s)
	%p = exp(-abs((x-mu)./s))./(2*s);
	%p = exp(-abs((x-mu)./s))./(2*s);
	p = laplacepdf(x,mu,s) + laplacepdf(x,-mu,s);
end


function f = gampdf_man(x,a,b)
	b = 1./b;
	f = exp(log(b).*a + log(x).*(a-1.0) - b.*x - gammaln(a));
end


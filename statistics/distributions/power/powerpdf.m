function pdf = powerpdf(x,a)
	 pdf = (a-1).*x.^(-a);
	 pdf(x < 1) = 0;
end


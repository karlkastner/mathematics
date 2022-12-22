function y=mybetainc(z,a,b,n)
	y = 0;
	for k=0:n
		y = y + z^a*(gamma(1-b+k)/gamma(1-b))/(factorial(k)*(a+k))*z^k;
	end
end

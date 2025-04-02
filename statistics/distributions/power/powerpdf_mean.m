function mu = powerpdf_mean(a,x0)
	mu = -(x0^(2 - a)*(a - 1))/(a - 2);
	%mu = (a-1)./(a-2);
	mu(a<2) = 0;
end


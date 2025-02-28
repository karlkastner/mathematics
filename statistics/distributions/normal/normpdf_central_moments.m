function m = normpdf_central_moments(sigma,k)
	m = sigma.^k.*double_factorial(k-1).*mod(k-1,2);
end


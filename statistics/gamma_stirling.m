function G = gamma_stirling(z)
	G = sqrt(2*pi).*z.^(z-1/2).*exp(-z);
end

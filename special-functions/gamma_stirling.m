% 2022-01-07 13:42:55.957389133 +0100
% Karl Kastner, Berlin
function G = gamma_stirling(z)
	G = sqrt(2*pi).*z.^(z-1/2).*exp(-z);
end

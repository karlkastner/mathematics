function [sk] = logn_skewness(lmu,lsd)
	sk = (exp(lsd.^2)+2).*sqrt(exp(lsd.^2)-1);
end


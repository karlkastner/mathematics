function R = cauchyacf(x,f0,s)
	R = cos(2*pi*f0*x).*exp(-(2*pi)*s*abs(x));
end


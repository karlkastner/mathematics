function xi = quantile_extrap(x,pi,dim)
	n = size(x,dim);
	x = sort(x,dim);
	p = (1:n)/(n+1);
	if (3==dim)
		xi = interp1(p,x,pi,'spline');
	end
end

function w = cauchypdf_width(x0,s,p0)
	if (nargin()<3)
		p0 = 0.5;
	end
	% syms x x0 s p0; solve(p0*cauchypdf(x0,x0,s) - cauchypdf(x,x0,s),x)
	w = 2*s*sqrt((1 - p0)/p0);
end


% 2024-06-29 11:49:32.647009759 +0200
function w = laplacepdf_width(x0,s,p0)
	if (nargin()<3)
		p0 = 0.5;
	end
	% syms x x0 s p0; solve(p0*laplacepdf(x0,x0,s) - cauchypdf(x,x0,s),x)
	w = -2*s*log(p0);
end


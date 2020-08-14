% Wed 13 May 21:17:39 +08 2020
function [Dx,Dy] = vanderd_2d(x,y,n)
	a    = length(x);
	b    = (n+1)*(n+2)/2;

	Dx = zeros(a,b);
	Dy = zeros(a,b);
	if (issym(x) || issym(y))
		Dx = sym(Dx);
		Dy = sym(Dy);
	end

	s = 0;
	m = 1;
	for idx=1:n
		s = s+idx;
		for jdx=1:idx
			Dx(:,s+jdx) = (idx+1-jdx).*x.^(idx-jdx).*y.^(jdx-1);
			Dy(:,s+jdx+1) = jdx*x.^(idx-jdx).*y.^(jdx-1);
		end
	end
end

% 2015-05-15 20:17:18.200830218 +0200
% Karl Kastner, Berlin

function mode = mode_man(x)
	fdx = isfinite(x);
	x = x(fdx);
%	n = length(x);
%	while (length(x) > 2)
%		l = min(x) + (max(x)-min(x))/3;
%		r = min(x) + 2/3*(max(x)-min(x));
%		ldx = x < l;
%		rdx = x > r;
%		if (sum(ldx) > sum(rdx))
%			x(rdx) = [];
%		else
%			x(ldx) = [];
%		end
%		[max(x) min(x)]
%		m = 0.5*(max(x)+min(x));
%		fdx = x > m;
%		n_ = sum(fdx);
%		if (n_ > n/2)
%			x = x(fdx);
%			n = n_;
%		else
%			x = x(~fdx);
%			n = n-n_;
%		end
%	end
%	if (n==2)
%	if (2 == length(x))
%		mode = mean(x);
%	else
%		mode = x;
%	end
	np = 5;
	n = length(x);
	m = round(sqrt(n));
%	h = 1/sqrt(n);
%	p = 0.5*h:h:1;
	p_ = cvec((0:m)/m);
	p  = cvec((1:m)/(m+1));
	q = quantile(x,p_);
	A = vander_1d(q,np);
	f = A \ p_
	A = vander_1d(q,np)
	g = flipud(cvec(polyder(flipud(f))));
	g2 = flipud(cvec(polyder(flipud(g))));
	roots(flipud(g2))
	plot(q,[A*f A(:,1:end-1)*g])
%	d = diff(q);
%	clf
%	plot(q(1:end-1),d)
%	[void mdx] = min(d);
%	mode = 0.5*(q(mdx)+q(mdx+1));
end


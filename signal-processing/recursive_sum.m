% 2015-06-22 17:22:56.585131862 +0200
% Karl Kastner, Berlin

function r = recsum(p,nmax)
	r = recsum_(p,0,nmax);
end

function r = recsum_(p,n,nmax)
	if (n <= nmax)
		r = p^n + p*recsum_(p,n+1,nmax) + recsum_(p,n+1,nmax);
		%r = p^n + p*recsum_(p,n+1,nmax);
		% 1/(1-p)
		%r = p^n + recsum_(p,n+1,nmax);
	else
	r = 0;
	end
end


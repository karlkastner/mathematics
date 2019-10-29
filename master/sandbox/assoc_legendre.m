% Thu Jan 19 02:01:29 MSK 2012
% Karl Kästner, Berlin
%
% 1/2^n 1/n! d^n/dx^n (x^2 - 1)^n
% test with matlab : mfun('P',n,x)
function f = assoc_legendre(l,m,x)
	f = (-1)^m * 1/2^l * 1/factorial(l) * (1 - x.^2).^(m/2) ...
		.*tail(x,0,l,l+m);
%k_max = floor(n/2);
%	f_    = 0.0;
%	for k=0:k_max
%		f_ = f_ + (-1)^k*(factorial(2*n - 2*k)*x.^(n-2*k))/(factorial(n-k)*factorial(n-2*k)*factorial(k)*2^n);
%	end
%	f_

	% d^k/dx^k ( x^m (x^2 - 1)^n )
	function f = tail(x,m,n,k)
		% end of recursion
		if (0 == k)
			f = x.^m.*(x.^2 - 1).^n;
		else
			f = m*tail(x,m-1,n,k-1) + 2*n*tail(x,m+1,n-1,k-1);
		end
	end
end


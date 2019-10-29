% Thu Jan 19 03:15:57 MSK 2012
% Karl Kästner
%
% L(n,k,x) = e^x x^-k 1/n! d^n/d^n ( e^-x x^(n+k) )

% todo k != 0 faulty
function [f_ f] = assoc_laguerre(n,k,x)
	%f = 1/factorial(n) * exp(x) .* x.^(-k) .* tail(x, n, n+k);
	f = 1/factorial(n) * x.^(-k) .* tail(x, n, n+k);
	% d^m/dx^m (x^n e^-x) = n d^m-1/dx^m-1 x^n-1 e^-x - d^m-1/dx^m-1 x^n e^-x
	function f = tail(x,n,m)
		if (0 == m)
			f = x.^n;
		else
			f = n*tail(x,n-1,m-1) - tail(x,n,m-1);
		end
	end

	f_ = 0;
	for m=0:n
		f_ = f_ + (-1)^m / (factorial(n-m)*factorial(k+m)*factorial(m)) * x.^m;
	end
	f_ = f_ * factorial(n+k);
end % function assoc_laguerre


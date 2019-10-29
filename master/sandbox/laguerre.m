% Dec 23 00:33 MSK 2011
% Karl Kästner, Berlin
%
% non-generalised laguerre polynomial
% L = 1/n! e^x d^n/dx^n (x^n exp(-x))
% does the same as mfun('L',x,n)
function f = laguerre(n,x)
	f = exp(x)./factorial(n).*tail(x,n,n).*exp(-x);
	% d^m/dx^m (x^n e^-x) = n d^m-1/dx^m-1 x^n-1 e^-x - d^m-1/dx^m-1 x^n e^-x
	function f = tail(x,n,m)
		if (0 == m)
			f = x.^n;
		else
			f = n*tail(x,n-1,m-1) - tail(x,n,m-1);
		end
	end
end % laguerre



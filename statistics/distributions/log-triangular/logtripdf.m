% Fri 23 Mar 12:15:10 CET 2018
%
%% probability density of the logarithmic triangular distribution
function f     = logtripdf(a,b,c,x)
	d = -b/2*log(a/c);

	% x < b
	f      = (log(x/a)./x)/(log(b/a)/b);
	fdx    = x > b;
	f(fdx) = (log(c./x(fdx))./x(fdx))/(log(c/b)/b);

	fdx    = (x<a) | (x>c);
	f(fdx) = 0;
	
	% normalize
	f      = f/d;
end


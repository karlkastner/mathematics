% Fri 23 Mar 12:15:10 CET 2018
%
%% density of the logarithmic triangular distribution
%
function f     = logtripdf(a,b,c,x)
	f      = (log(x) - log(a))/(log(b)-log(a));
	fdx    = x<a;
	f(fdx) = 0;
	fdx    = x>b;
	f(fdx) = (log(c) - log(x(fdx)))/(log(c)-log(b));
	f(x>c) = 0;
	%den    = (b - c - b*log(b) + b*log(c))/((log(b) - log(c)));
	den =     (b*(log(c) - log(b) + 1))/((log(b) - log(c))) ...
		- (a - b - b*log(a) + b*log(b))/((log(a) - log(b))) ...
	      - c/((log(b) - log(c)));
	f      = f/den;
end


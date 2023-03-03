% Fri 23 Mar 12:20:14 CET 2018
%
%% cumulative distribution of the logarithmic triangular distribution
%
function F = logtricdf(a,b,c,x)

	d = -b/2*log(a/c);

	% x < b
	F      = (b*(log(a) - log(x)).^2)/(2*d*(log(b) - log(a)));
	fdx    = x>b;
 	Fb     = (b*(log(c) - log(b)))/(2*d) - (b*(log(a) - log(b)))/(2*d);
	F(fdx) = Fb + (b*(log(c) - log(x(fdx))).^2)/(2*d*(log(b) - log(c)));

	fdx    = x<a;
	F(fdx) = 0;
	fdx    = x>c;
	F(fdx) = 1;
end % logtripdf


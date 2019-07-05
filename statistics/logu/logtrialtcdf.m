% Fri 23 Mar 12:20:14 CET 2018
%
%% pdf of a logarithmic triangular distribution
%
function F = logtripdf(a,b,c,x)
	F = -(a - x - x*log(a) + x.*log(x))/((log(a) - log(b)));
	fdx = x<a;
	F(fdx) = 0;
	fdx = x>b;
	%F(fdx) = -(a - x(fdx) - x(fdx)*log(a) + x(fdx).*log(x(fdx)))/((log(a) - log(b)));
	F(fdx) =   (b*(log(c) - log(b) + 1))/((log(b) - log(c))) ...
		 - (a - b - b*log(a) + b*log(b))/((log(a) - log(b))) ...
                 - (x(fdx).*(log(c) - log(x(fdx)) + 1))/((log(b) - log(c)));
	%den = (b - c - b*log(b) + b*log(c))/((log(b) - log(c)));
	den    =   (b*(log(c) - log(b) + 1))/((log(b) - log(c))) ...
                 - (a - b - b*log(a) + b*log(b))/((log(a) - log(b))) ...
                 - c/((log(b) - log(c)));
	F = F/den;
	fdx = x>c;
	F(fdx) = 1;
end % logtripdf


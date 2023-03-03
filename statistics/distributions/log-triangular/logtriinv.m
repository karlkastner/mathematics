% Fri 23 Mar 13:35:25 CET 2018
%
%% invere of the logarithmic triangular distribution
function x = logtriinv(a,b,c,F)
	den    = -(b*(log(a) - log(c)))/(2);

	Fb     = -(b*(log(a) - log(b)))/(2*den);
	x      = a*exp((2^(1/2)*F.^(1/2)*den^(1/2)*(log(b) - log(a))^(1/2))/b^(1/2));
	fdx    = F > Fb;
	x(fdx) = c*exp(-((log(b) - log(c))*(2*F(fdx)*den + b*log(a) - b*log(c))).^(1/2)/b^(1/2));

	x(F<0) = NaN;
	x(F>1) = NaN;
end


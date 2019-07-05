% 2015-02-11 14:24:37.740894622 +0100
% Karl Kastner, Berlin
%
%% correction factor for standard error of the mean of n ar1-correlated iid samples
%
% function [f f_simple] = f_correlation(rho, n)
function [f f_simple] = f_correlation(rho, n)
	d         = rho*((n-1) - n*rho + rho.^n)/(1-rho)^2;
	%d        = ((n-1)*rho - n*rho*rho +rho.^(n+1))/(1-rho)^2;
	f2        = (n-1)*(n + 2*d)./(n*(n-1)-2*d);
	%f2       = (n + 2*d)./(n-2*d./(n-1));
	%f2       = (1 + 2*d./n)./(1-2*d./(n.*(n-1)));
	f         = sqrt(f2);

	% simplified correction factor, limit case for n->infty
	f_simple2 = (1+rho)./(1-rho);
	f_simple  = sqrt(f_simple2);
end


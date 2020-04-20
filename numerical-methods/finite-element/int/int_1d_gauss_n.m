% Mon 30 Mar 10:41:19 +08 2020
function [wi,xi] = int_1d_gauss_n(n)
	% coefficients of the legendre polynomial
	x0    = linspace(-1,1,n+1)';
	A     = fliplr(vander_1d(x0,n));
	y0    = legendreP(n,x0);
	c_Pn  = A \ y0;
	% derivative of the legendre Polynomial
	c_dPn = polyder(c_Pn);

	xi = roots(c_Pn);
	wi = 1./((1-xi.^2).*polyval(c_dPn,xi).^2);
	% transform to 0, 1
	xi = 0.5*(1+xi);
	% barycentric coordinates	
	xi = [xi,1-xi];
end

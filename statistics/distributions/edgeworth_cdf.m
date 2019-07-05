% Do 11. Feb 18:11:02 CET 2016
%
%% edgeworth expansion of an unknown cumulative distribution
%% with mean mu, standard deviation sigma, and third and fourth cumulants
%% c.f. Rao 2010
function F = edgeworth_cdf(mu,sigma,c3,c4,x)
	% standardise
	x  = (x-mu)/sigma;
	p1 = -1/6*c3*(x.^2-1);
	% these are hermite polynomials
	p2 = -x.*(1/24*c4*(x.^2-3) ...
	     + 1/72*c3^2*(x.^4 - 10*x.^2 + 15));
	F  = normcdf(x) + normpdf(x).*(p1 + p2);
end


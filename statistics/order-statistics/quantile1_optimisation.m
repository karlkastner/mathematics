% 2016-03-01 13:57:53.913565257 +0100
% Karl Kastner, Berlin

function q = quantile1_optimisation(x,p)
	q = mean(x);
	%q = lsqnonlin(@(x) objective(x,q,p),q);
	q = fminsearch(
	[h g f] = hessian(@(x) objective(x,q,p),q)
end
function f = objective(x,q,p)
	fdx = (x < q);
	u = 1;
	rho = u*(p - 1.0*fdx);
	% f = (1-p)*sum(abs(x(fdx)-q)) + p*sum(abs(x(~fdx)-q));
	f = rho.*(x - q);
end

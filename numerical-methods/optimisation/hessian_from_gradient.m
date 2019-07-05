% Thu  8 Mar 11:43:50 CET 2018
%
%% numerical hessian from gradient
function [h] = hessian_from_gradient(gfun,x0)
	p = 1/4;
	dx = max(abs(x0).*eps.^p,eps.^p);
	h = grad(gfun,x0,dx);
%	h2 = grad(gfun,x0);
end


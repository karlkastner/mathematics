% 2024-05-24 17:06:08.578822484 +0200
function [mu,gamma] = cauchycdf_fit(x,c)
	x = cvec(x);
	c = cvec(c);
	nx = length(x);
	A = [ones(nx,1),x];
	% c = 0.5+1/pi*atan((x-mu)/gamma);
	rhs = tan(pi*(c-0.5));
	par = A \ rhs;
	% tan(pi*(c-0.5)) = x/g - mu/g
	gamma = 1./par(2);
	mu = -par(1)/par(2);
end


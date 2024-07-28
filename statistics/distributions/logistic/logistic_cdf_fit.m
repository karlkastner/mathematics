function [mu,a] = logistic_cdf_fit(x,c)
	x = cvec(x);
	c = cvec(c);
	% log(-1./(c-1)-1) = a*x - a*mu
	par = [ones(length(x),1),x]\log(-1./(c-1)-1);
	mu = -par(1)/par(2);
	a  = par(2);
end

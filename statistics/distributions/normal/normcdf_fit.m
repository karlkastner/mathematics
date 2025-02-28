function [mu,sd] = normcdf_fit(x,c)
	c = cvec(c);
	% c = normcdf(x,a,b)
	% norminv(c) = (x-a)/b
	% a + b*norminv(v) = x
	par = [ones(length(c),1),norminv(c)] \ cvec(x);
	mu = par(1);
	sd = par(2)
end

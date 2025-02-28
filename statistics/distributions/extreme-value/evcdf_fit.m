function [mu,b] = gumbel_cdf_fit(x,F)
	% F = exp(-exp(-(x-mu)/b))
	% -log(-log(F)) = x/b - mu/b 
	x = cvec(x);
	F = cvec(F);
	n = length(x);
	% alternative matlab definition: (as in wolfram)
	%  F = 1-exp(-exp((x-mu)/b));
	A  = [ones(n,1),x];
	c  = A \ log(-log(1-F));
	b  = 1./c(2);
	mu = -c(1)/c(2);
end


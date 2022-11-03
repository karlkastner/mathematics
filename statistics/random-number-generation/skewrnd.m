% 2015-08-03 14:09:44.313605435 +0200
% Karl Kastner, Berlin
%
%% random numbers of the skew normal distribution
% function [x mu sd sk] = skewnrnd(mu,sd,sk,n)
function [x, mu, sd, sk] = skewnrnd(mu,sd,sk,n,n2)
	a = skewness2param(sk);
	if (nargin()>4)
		n(2) = n2;
	end
	if (length(n)==1)
		n(2) = 1;
	end
	% sample
%	% skew normal - this does not work, use algorithm
	% x = 2*x.*normcdf(s*x);
	% X = randn(n,m);
	% X = X.*normcdf(a*X);
	z = zeros(n);
	m = prod(n);
	id=(1:m);
	while (m > 0)
		x   = randn(m,1);
		y   = randn(m,1);
		sdx = x < a*y;
		z(id(sdx)) = y(sdx);
		id(sdx)    = [];
		m          = length(id);
	end
	x = z;
	% "normalisation" (except for skewneww)
	[sk, sd_, mu_] = skewpdf_central_moments(a);
	x = (x-mu_)*(sd/sd_)+mu;
	%x = x*sd/sd_-(mu_mu_);
end



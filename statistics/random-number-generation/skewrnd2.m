% Do 11. Feb 18:53:50 CET 2016
% Karl Kastner, Berlin
%
%% random numbers of the skew normal distribution
function X = skewrnd(a,n,m)
	if (nargin() < 3)
		m = 1;
	end
	X = randn(n,m);
	X = X.*normcdf(a*X);
end

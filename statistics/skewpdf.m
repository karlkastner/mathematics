% 2017-02-27 10:13:18.326037478 +0100
% Karl Kastner, Berlin
%
%% skew-normal distribution
%% c.f. Azzalini 1985
%function y = skewpdf(x,mu,sd,sk)
function y = skewpdf(x,mu,sd,sk)

	a = skewness2param(sk);
	[sk_, sd_, mu_] = skewparam_to_central_moments(a);
	
	% translate and scale
	x = ((x-mu)*sd_/sd+mu_);
	
	y = 2*sd_/sd*normcdf(a*x).*normpdf(x);
end % skewpdf


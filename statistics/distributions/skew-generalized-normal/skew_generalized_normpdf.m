% Mon 23 Nov 12:55:32 +08 2020
% Karl Kastner, Berlin
%% c.f. Arellano-Valle 2004
%
% note : mu,sd,sk and ku are parameters here,
%        actual mean,sd sk and ku will differ,
%        exact transformation is not possible, due to lack of exact solution
%        for odd moments
function y = skew_generalized_normpdf(x,mu,sd,l1,l2)

if (0)
	a  = skewness2param(sk);
	l1 = a;
	l2 = ku;
	if (l2<0)
		error('l2 must be larger 0');	
	end
	[sk_, sd_, mu_] = skewparam_to_central_moments(a);

	% translate and scale
	x = ((x-mu)*sd_/sd+mu_);

	y = 2*sd_/sd*normcdf(a*x).*normpdf(x);
else
	z  = (x-mu)./sd;

%	l1 = sk;
%	l2 = ku;
%	l1 = tan(l1*(0.5*pi));
%	l1 = l1/(1-abs(l1))

	% note : this can be reparametrized as l2' = l1*l2,
	% so that the skewness and kurtosis can be computed with a
	% single parameter l2', given that the parameters of the skew-normal distribution are used
	y  = 2/sd*normcdf(l1*z./sqrt(1+l2*z.^2)).*normpdf(z);
	%y   = 2/sd*normcdf(l1*z./sqrt(1+(l1*l2*z).^2)).*normpdf(z); this was active
	%y  = 2/sd*normcdf(l1*z./sqrt(1+(l2*z).^2)).*normpdf(z);
end		

end % skew_generalized_normpdf


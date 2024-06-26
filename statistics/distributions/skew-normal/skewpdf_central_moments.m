% 2015-08-03 14:09:44.313605435 +0200
% Karl Kastner, Berlin
function [sk, sd, mu] = skewpdf_central_moments(a)
	d = a/sqrt(1+a^2);
	b = sqrt(2/pi);
	% mean
	mu = b*d;
	% standard deviation
	sd = sqrt((1-(b*d)^2));
	% skewness
	sk = 0.5*(4-pi)*mu^3/sd^3;
end



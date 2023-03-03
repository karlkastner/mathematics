% 2017-02-27 10:13:18.326037478 +0100
% Karl Kastner, Berlin

function [sk, sd, mu] = skewparam_to_central_moments(a)
	d = a/sqrt(1+a^2);
	b = sqrt(2/pi);
	% mean
	mu = b*d;
	% standard deviation
	sd = sqrt((1-(b*d)^2));
	% skewness
	sk = 0.5*(4-pi)*mu^3/sd^3;
end


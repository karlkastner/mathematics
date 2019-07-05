% 2017-02-27 10:13:18.326037478 +0100
% Karl Kastner, Berlin
%
%% skew-normal distribution
%% c.f. Azzalini 1985
%function y = skewpdf(x,mu,sd,sk)
function y = skewpdf(x,mu,sd,sk)
%	x = (x-mu)/sd;

	a = param(sk);
	[sk_ sd_ mu_] = central_moments(a);
	
	x = ((x-mu)*sd_/sd+mu_);
	
	y = 2*sd_/sd*normcdf(a*x).*normpdf(x);
	y = y; % sd?
end

function a = param(sk)
	if (abs(sk) > 1-sqrt(eps))
		error('skewness has to be strictly between -1 and 1');
	end
	sk23 = abs(sk).^(2/3);
	% solve for delta
	delta = sqrt(sk23/(sk23*2/pi + (0.5*(4-pi))^(2/3)*2/pi));
	a = sign(sk)*delta/sqrt(1-delta^2);
end


function [sk sd mu] = central_moments(a)
	d = a/sqrt(1+a^2);
	b = sqrt(2/pi);
	% mean
	mu = b*d;
	% standard deviation
	sd = sqrt((1-(b*d)^2));
	% skewness
	sk = 0.5*(4-pi)*mu^3/sd^3;
end


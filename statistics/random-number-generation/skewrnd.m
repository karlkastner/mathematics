% 2015-08-03 14:09:44.313605435 +0200
% Karl Kastner, Berlin
%
%% random numbers of the skew normal distribution
% function [x mu sd sk] = skewnrnd(mu,sd,sk,n)
function [x, mu, sd, sk] = skewnrnd(mu,sd,sk,n)
	a = param(sk);

	% sample
%	% skew normal - this does not work, use algorithm
	% x = 2*x.*normcdf(s*x);	fdx = true(n,1);
	fdx=(1:n);
	m = n;
	while (m > 0)
		x = randn(m,1);
		y = randn(m,1);
		sdx = x < a*y;
		z(fdx(sdx)) = y(sdx);
		fdx(sdx) = [];
		m = length(fdx);
	end
	x = z;
	% "normalisation" (except for skewneww)
	[sk sd_ mu_] = central_moments(a);
	x = (x-mu_)*(sd/sd_)+mu;
	%x = x*sd/sd_-(mu_mu_);
end

function a = param(sk)
	if (abs(sk) > 1-sqrt(eps))
		error('skewness has to be between -1 and 1');
	end
	sk23 = abs(sk).^(2/3);
	% solve for delta
	delta = sqrt(sk23/(sk23*2/pi + (0.5*(4-pi))^(2/3)*2/pi));
%	delta = sqrt(0.5*pi*sk23/(sk23+(0.5*(4-pi))^(2/3)))
	a = sign(sk)*delta/sqrt(1-delta^2);
%	a = fzero(@(a) sk - central_moments(a), 0);
%	d = a/sqrt(1+a^2)
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


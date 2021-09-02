% Sat 28 Aug 14:15:01 CEST 2021
% wrapped normal distribution to the unit circle
% stephens
function f = wnormpdf(theta,mu,sigma,m)
	f = 0;
	for k=-m:m
		f = f + exp(-(theta - mu + 2*pi*k).^2./(2*sigma.^2));
	end
	f = 1./(sigma*sqrt(2*pi)).*f;
end


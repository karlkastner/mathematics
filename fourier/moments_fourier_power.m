% 2022-04-09 14:42:42.274759731 +0200
% moments of (1 + cos(2*pi*x/L))^p
function m = moments_fourier_power(p,L)
	if (nargin()<2)
		L = 1;
	end
	m = L*(gamma(2*p+1)./(gamma(p+1).^2))./2.^p;	
end


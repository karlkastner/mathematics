% 2022-04-09 14:42:42.274759731 +0200
% variance of (1 + cos(2*pi*x/L))^p
function s2 = var_fourier_power(p,varargin)
	m1 = moments_fourier_power(p,varargin{:});
	m2 = moments_fourier_power(2*p,varargin{:});
	s2 = m2 - m1.^2;
end

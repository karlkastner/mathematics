% 2022-04-09 14:42:42.274759731 +0200
% variance of (1 + cos(2*pi*x/L))^p
function s = std_fourier_power(varargin)
	s2 = var_fourier_power(varargin{:});
	s  = sqrt(s2);
end


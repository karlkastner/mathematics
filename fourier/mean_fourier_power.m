% 2022-04-09 14:42:42.274759731 +0200
% mean of (1 + cos(2*pi*x/L))^p
function m = mean_fourier_power(p,varargin)
	m = moments_fourier_power(p,varargin{:});
end


% 2022-04-03 10:13:22.054340576 +0200
function w = caesaro_weight(varargin)
	[f, T, mask, N] = fourier_axis(varargin{:});
	% TODO n or n+1 here?
	n = length(N);
	w = (1-2*abs(N)/(n+1)).^(1);
end

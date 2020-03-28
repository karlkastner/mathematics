% 2016-02-29 15:14:30.485135195 +0100
% Karl Kastner, Berlin
function h = plotcdf(x,m,varargin)
	x = sort(x(isfinite(x)));
	n = length(x);
	p = (1:n)/(n+1);
	% TODO resample
	h = plot(x,p,varargin{:});
end


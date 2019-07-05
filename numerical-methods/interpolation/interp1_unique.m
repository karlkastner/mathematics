% Wed 17 Aug 18:05:18 CEST 2016
% Karl Kastner, Berlin
%% matlab fails to interpolate, when x values are not unique
%% this function makes the values unique before use
function yi = interp1_unique(x,y,xi,varargin)
	% exclude NaN
	fdx        = isfinite(x) & isfinite(y);
	x          = x(fdx);
	y          = y(fdx);
	if (length(x)>1)
		% make x unique and sort
		[xu id jd] = unique(x);
		% workaround matlab bug
		x0 = xu(1);
		yi = interp1(xu-x0,y(id),xi-x0,varargin{:});
	else
		yi = NaN(size(xi));
	end
end


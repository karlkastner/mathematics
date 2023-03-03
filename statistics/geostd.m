% Thu 26 Mar 12:49:04 +08 2020
% function sd = geostd(x)
% c.f. wikipedia
function gsd = geostd(x,varargin)
	n   = length(x);
	% sg = exp(sqrt( sum( (log(x) - mean(log(x))).^2 )./(n-1)))
	gsd = exp(std(log(x),varargin{:}));
	%sd  = geomean(x).*std(log(x));
end


% Fr 15. Mai 14:25:31 CEST 2015
% Karl Kastner, Berlin
%
%% norm of a vector, skips nan-values
function n = nannorm(x,varargin)
	fdx = isfinite(x);
	n = norm(x(fdx),varargin{:});
end


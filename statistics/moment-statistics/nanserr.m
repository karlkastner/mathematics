% Thu Nov 27 22:08:26 CET 2014
% Karl Kastner, Berlin 
%
%% standard error of x with respect to mean when x contains nan values
function s = nanserr(x,dim)
	if (isempty(x))
		s = NaN(class(x));
	elseif (nargin < 2 || isempty(dim))
		%dim = 1;
		s = nanstd(x)./sqrt(sum(isfinite(x)));
	else
		s = nanstd(x,[],dim)./sqrt(sum(isfinite(x),dim));
	end
end

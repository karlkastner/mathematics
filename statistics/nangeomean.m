% 2017-06-28 12:57:43.177797184 +0200
%
% geometric mean excluding nan-values
function x = nangeomean(x,varargin)
	x = exp(nanmean(log(x),varargin{:}));
%	if (nargin()<2)
%		dim = 1;
%		if (isvector(x))
%			x = cvec(x);
%		end
%	end
%	n = size(x,dim);
%	s = size(x);
%	s(dim) = 1;
%	y = NaN(s);
%	for idx=1:n
%		if (1==dim)
%			y(1,:) = 
%	end
end


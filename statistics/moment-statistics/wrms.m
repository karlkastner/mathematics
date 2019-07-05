% Wed  9 Aug 11:56:50 CEST 2017
%
%% weighted root mean square error
% varargin can be dimension
%
function wrms = wrms(w,x,varargin)
	if (isvector(x))
		x = cvec(x);
	end
	if (isvector(w))
		w = cvec(w);
	end
	wrms = sqrt(sum(bsxfun(@times,w.^2,x.^2),varargin{:})./sum(w.^2,varargin{:}));
end

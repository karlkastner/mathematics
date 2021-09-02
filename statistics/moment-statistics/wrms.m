% Wed  9 Aug 11:56:50 CEST 2017
%
%% weighted root mean square
% varargin can be dimension
%
%     mu = 1/(sum w) sum w x
%      e = x - mu
%
% function wrms = wrms(w,x,varargin)
function wrms = wrms(w,x,varargin)
	if (isvector(x))
		x = cvec(x);
	end
	if (isvector(w))
		w = cvec(w);
	end
	wms  = sum(bsxfun(@times,w,x.^2),varargin{:})./sum(w,varargin{:});
	wrms = sqrt(wms);
end

% 2015-04-07 12:26:07.997739673 +0200
% Karl Kastner, Berlin
%
%% weighted mean
%%
%% min_x sum w (x-mu)^2 => mu = sum(wx)/sum(w)
%%
%% varargin can be dim
%% function [mu serr] = wmean(w,x)
%
% TODO should negative weights be interpretred as sign change,
% so that sw = sum(|w|) ? -> no, think of sample being added twice
function [mu, serr] = wmean(w,x,varargin)
	if (isvector(x) && nargin() < 3)
		x = cvec(x);
	end
	n2 = size(x,2);
	if (isvector(w))
		w = repmat(cvec(w),1,n2);
	end

	% scale values by weights
	wx       = w.*x;

	% number of samples after weighing (does not equal effective sample size)
	sw  = sum(w,varargin{:});

	% weighted mean
	mu  = sum(wx,varargin{:})./sw;

	if (nargout() > 1)
		wdx2  = w.*bsxfun(@minus,x,mu).^2;
		% number of degree of freedoms (this is the effective sample size)
		ni    = sw.^2./sum(w.^2);
		% error variance
		s2    = 1./(sw.*(ni-1)).*sum(wdx2);
		% standard error
		serr = sqrt(s2);
	end % if
end % wmean


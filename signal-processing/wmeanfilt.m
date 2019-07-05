% Di 2. Feb 10:11:54 CET 2016
% Karl Kastner, Berlin
%% mean filter with window
% function [me_ serr_] = wmeanfilt(w,x,detrend,varargin)
function [me_ serr_] = wmeanfilt(w,x,detrend,varargin)
	if (nargin() < 3)
		detrend = false;
	end

	if (isvector(x))
		x=cvec(x);
	end

	nf  = length(w);

	% make even length filter windows uneven
	% if this is not done, the output series is shifted by 1
	if (0 == mod(nf,2))
		w(end+1) = 0;
		w(2:end) = w(2:end)+w(1:end-1);
	end

	% normalise filter weights
	w = w/sum(w);

	nx    = size(x);

	% allocate memory
	me   = zeros(nx);
	serr = zeros(nx);
	for idx_=1:size(x,2);
		[me_(:,idx_) serr_(:,idx_)] = wmeanfilt_(w,x(:,idx_),detrend,varargin{:});
	end % for
end

function [me serr] = wmeanfilt_(w,x,detrend,varargin)

	nx  = length(x);
	nf  = length(w);
	nfl = floor(nf/2);
	nfr = nf-nfl-1;

	% allocate memory
	me   = NaN(nx,1);
	serr = NaN(nx,1);

	% left end
	for idx=1:nfl
		r = nfr+idx;
		jdx=(1:r)';
		w_ = w(nfl-idx+2:nf); % 1:nf-1); ?
		fdx = (isfinite(x(jdx)));
		w_ = w_(fdx);
		w_ = w_/sum(w_);
		jdx = jdx(fdx);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = PolyOLS.detrend(jdx-idx,x(jdx),w_);
		end
		if (~isempty(x_))
		[me(idx) serr(idx)] = wmean(w_,x_,varargin{:});
		end
	end % for idx

	% centre
	% TODO, this can be implemented faster with convolution
	parfor idx=nfl+1:nx-nfr
		%l = max(idx-nfl,1);
		%r = min(n,idx+nfr);
		l = idx-nfl;
		r = idx+nfr;
		jdx = (l:r)';
		w_  = w;
		fdx = (isfinite(x(jdx)));
		w_  = w_(fdx);
		w_  = w_/sum(w_);
		jdx = jdx(fdx);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = PolyOLS.detrend(jdx-idx,x(jdx),w_);
		end
		if (~isempty(x_))
		[me(idx) serr(idx)] = wmean(w_,x_,varargin{:});
		end
	end % for idx

	% right end
	for idx=nx-nfr:nx
		l = idx-nfl;
		jdx=(l:nx)';
		w_ = w(1:nfl+nx-idx+1);
		fdx = (isfinite(x(jdx)));
		w_ = w_(fdx);
		w_ = w_/sum(w_);
		jdx = jdx(fdx);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = PolyOLS.detrend(jdx-idx,x(jdx),w_);
		end
		if (~isempty(x_))
		[me(idx) serr(idx)] = wmean(w_,x_,varargin{:});
		end
	end % for idx
end % wmeanfilt_



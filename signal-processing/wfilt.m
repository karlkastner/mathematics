% Di 2. Feb 10:11:54 CET 2016
% Karl Kastner, Berlin
%% filter with window
function varargout = wfilt(@func,w,x,varargin)
	if (nargin < 3)
		P = [];
	end

	n    = length(x);

	% allocate memory
	me   = NaN(n,1);
	me_s = NaN(n,1);
	me_l = NaN(n,1);
	me_u = NaN(n,1);

	nf  = length(w);
	nfl = floor(nf/2);
	nfr = nf-nfl-1;
	if (n >= nf)
	% left end
	for idx=1:nfl
		l = 1;
		r = nfr+idx;
		jdx=1:r;
		w_ = w(nfl-idx+1:nf);
		w_ = w_/sum(w_);
		[me(idx) me_s(idx) me_l(idx) me_u(idx)] = func(w_,x(jdx),varargin{:});
	end
	% centre
	parfor idx=nfl+1:n-nfr
		%l = max(idx-nfl,1);
		%r = min(n,idx+nfr);
		l = idx-nfl;
		r = idx+nfr;
		jdx = l:r;
		%w = ones(length(jdx),1);
		%(idx,:) = wmedian(w,x(jdx),p);
		[me(idx) me_s(idx) me_l(idx) me_u(idx)] = func(w,x(jdx),varargin{:});
	end
	% right end
	for idx=n-nfr:n
		l = idx-nfl;
		r = n;
		jdx=l:r;
		w_ = w(1:nfl+n-idx+1);
		w_ = w_/sum(w_);
		%y(idx,:) = wmedian(w_,x(jdx),p);
		[me(idx) me_s(idx) me_l(idx) me_u(idx)] = func(w_,x(jdx),varargin{:});
	end
	end
end


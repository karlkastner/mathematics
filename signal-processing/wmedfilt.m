% Di 2. Feb 10:11:54 CET 2016
% Karl Kastner, Berlin
%% median filter with window
function [me me_s me_l me_u] = wmedfilt(w,x,detrend,varargin)
	if (nargin()<3)
		detrend = false;
	end

	if (isvector(x))
		x=cvec(x);
	end

	theil = Theil();
	nf  = length(w);
		
	% make even length filter windows uneven
	% if this is not done, the output series is shifted by 1
	if (0 == mod(nf,2))
		w(end+1) = 0;
		w(2:end) = w(2:end)+w(1:end-1);
		nf = nf+1;
		% length(w);
		%w = interp1((0:nf-1)/(nf-1),w,(0:nf)/nf,'cubic');
	end

	% limit filter length to data length
	% TODO, this conflicts with even / uneven
%	if (nf > nx)
%		w = interp1((1:nf)'/(nf+1),w,(1:nx)/(nx+1));
%		nf = nx;
%	end
%	if (nf > nx)
%		% clip
%		d  = nf-n;
%		d  = ceil(d/2);
%		w  = w(d+1:end-d);
%		nf = length(w);
%	end


	nx    = size(x);

	% allocate memory
	me   = NaN(nx);
	me_s = NaN(nx);
	me_l = NaN(nx);
	me_u = NaN(nx);
	for idx=1:size(x,2);
	        [me(:,idx) me_s(:,idx) me_l(:,idx) me_u(:,idx)] = wmedfilt_(w,x(:,idx),detrend,varargin);
	end % for

end

function [me me_s me_l me_u] = wmedfilt_(w,x,detrend,varargin)
	fun = @(varargin) nanwmedian(varargin{:});
	fun = @(varargin) wmedian(varargin{:});

	theil = Theil();

	nx = length(x);
	nf = length(w);	

	% normalise filter weights
	w = w/sum(w);

	nfl = floor(nf/2);
	nfr = nf-nfl-1;

	% left end
	for idx=1:nfl
	%parfor idx=1:nfl
		r   = nfr+idx;
		jdx = (1:r)';
		w_  = w(nfl-idx+2:nf);
		w_  = w_/sum(w_);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = theil.detrend(jdx-idx,x(jdx),w_);
		end
		[me(idx), me_s(idx), me_l(idx), me_u(idx)] = fun(w_,x_,varargin{:});
	end % for idx

	% centre
	for idx=nfl+1:nx-nfr
	%parfor idx=nfl+1:nx-nfr
		l = idx-nfl;
		r = idx+nfr;
		jdx = (l:r)';
		% do only detrend interior values when flag is set
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = theil.detrend(jdx-idx,x(jdx),w);
		end
		[me(idx), me_s(idx), me_l(idx), me_u(idx)] = fun(w,x_,varargin{:});
	end % for idx

	% right end
	for idx=nx-nfr:nx
	%parfor idx=nx-nfr:nx
		l = idx-nfl;
		jdx=(l:nx)';
		w_ = w(1:nfl+nx-idx+1);
		w_ = w_/sum(w_);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = theil.detrend(jdx-idx,x(jdx),w_);
		end
		[me(idx), me_s(idx), me_l(idx), me_u(idx)] = fun(w_,x_,varargin{:});
	end % for idx
end % wmedfilt


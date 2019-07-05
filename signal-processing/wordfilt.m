% Di 2. Feb 10:11:54 CET 2016
% Karl Kastner, Berlin
%% weighted order filter
function y = wordfilt(w,x,p,detrend)
	if (nargin() < 4)
		detrend = false;
	end
	if (isvector(x))
		x=cvec(x);
	end
	y = zeros(size(x,1),size(x,2),length(p));
	for idx =1:size(x,2)
		y(:,idx,:) = wordfilt_(x(:,idx));
	end
	if (size(x,2) == 1)
		y = squeeze(y);
	end

function y = wordfilt_(x)

	nx = length(x);
	nf  = length(w);

	% allocate memory
	y = NaN(nx,length(p));

	% make even length filter windows uneven
	if (0 == mod(nf,2))
		w(end+1) = 0;
		w(2:end) = w(2:end)+w(1:end-1);
		nf = length(w);
	end

	% limit filter length to data length
	% TODO, this conflicts with even / uneven
	if (nf > nx)
		w = interp1((1:nf)'/(nf+1),w,(1:nx)/(nx+1));
		nf = nx;
	end

	% normalise filter weights
	w = w/sum(w);

	nfl = floor(nf/2);
	nfr = nf-nfl-1;

	% left end
	parfor idx=1:nfl
		r  = nfr+idx;
		jdx=(1:r)';
		w_ = w(nfl-idx+2:nf); % 1:nf ?
		w_ = w_/sum(w_);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = Theil.detrend(jdx-idx,x(jdx),w_);
		end
		y(idx,:) = wquantile(w_,x_,p);
	end % for idx

	% centre
	parfor idx=nfl+1:nx-nfr
		l = idx-nfl;
		r = idx+nfr;
		jdx = (l:r)';
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = Theil.detrend(jdx-idx,x(jdx),w);
		end
		y(idx,:) = wquantile(w,x_,p);
		%y(idx,:) = wquantile(w,x(jdx),p);
	end % for idx

	% right end
	parfor idx=nx-nfr:nx
		l = idx-nfl;
		jdx=(l:nx)';
		w_ = w(1:nfl+nx-idx+1);
		w_ = w_/sum(w_);
		if (0 == detrend)
			x_ = x(jdx);
		else
			x_ = Theil.detrend(jdx-idx,x(jdx),w_);
		end
		y(idx,:) = wquantile(w_,x_,p);
	end % for idx
end % wordfilt
end

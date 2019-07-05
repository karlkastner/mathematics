% Do 11. Feb 19:37:18 CET 2016
% Karl Kastner, Berlin
%% weighed order filter
function y = wordfilt(w,x,p)
	n = length(x);
	y = NaN(n,length(p));
	nf  = length(w);
	nfl = floor(nf/2);
	nfr = nf-nfl-1;
	if (n >= nf)
	% left end
	for idx=1:nfl
		l = 1;
		r = nfr+idx;
		jdx=1:r;
		w_ = w(nfl-idx+1:nf-1);
		w_ = w_/sum(w_);
		x_ = x(jdx);

		mu = wmean(w_,x_);
		sd = wstd(w_,x_);
		sk = wskew(w_,x_);
		ku = wkurt(w_,x_);
		y(idx,:) = edgeworth_quantile(mu,sd,[],[],p,sk,ku,1);
	end
	% centre
	for idx=nfl+1:n-nfr
		%l = max(idx-nfl,1);
		%r = min(n,idx+nfr);
		l = idx-nfl;
		r = idx+nfr;
		jdx = l:r;
		%w = ones(length(jdx),1);
		x_ = x(jdx);
		w_ = w/sum(w);

		mu = wmean(w_,x_);
		sd = wstd(w_,x_);
		sk = wskew(w_,x_);
		ku = wkurt(w_,x_);
		y(idx,:) = edgeworth_quantile(mu,sd,[],[],p,sk,ku,1);
	end
	% right end
	for idx=n-nfr:n
		l = idx-nfl;
		r = n;
		jdx=l:r;
		w_ = w(1:nfl+n-idx+1);
		w_ = w_/sum(w_);
		x_ = x(jdx);

		mu = wmean(w_,x_);
		sd = wstd(w_,x_);
		sk = wskew(w_,x_);
		ku = wkurt(w_,x_);
		y(idx,:) = edgeworth_quantile(mu,sd,[],[],p,sk,ku,1);
	end
	end
end


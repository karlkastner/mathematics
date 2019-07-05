% Di 7. Apr 16:32:28 CEST 2015
% Karl Kastner, Berlin
%
%% weighted quantile, skips nan values 
% TODO negative weights are not too uninteresting
function [q] = nanwquantile(w,x,p,varargin)
	% cubic
	method = 'linear';
	nx     = size(x,2);
	np     = length(p);
	q      = NaN(np,nx);

	if (isvector(w))
		w = repmat(cvec(w),1,nx);
	end
	flag   = isnan(x) | isnan(w);
	w(flag) = 0;
	x(flag) = inf;

	% sort x and weigths
	[xs sdx] = sort(x,1);
	for idx=1:nx
		ws(:,idx) = w(sdx(:,idx),idx);
	end
	% integrate weights
	% TODO, this is non-central
	sws = sum(w);
	ws = bsxfun(@times,cumsum(ws),1./sws);

	q = NaN(length(p),size(x,2));
	for idx=1:nx
		if (sws(idx) > 1) %nf(idx)>1 && size(x,1)>1)
			% exclude duplicates
			[ws_ udx] = unique(ws(:,idx));
			q(:,idx)  = interp1(ws_,xs(udx,idx),p,method,NaN);
		end
	end
end % nanwquantile


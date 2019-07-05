% Di 7. Apr 16:32:28 CEST 2015
% Karl Kastner, Berlin
%
%% weighted quantile
% assume finite weights
function [q] = wquantile(w,x,p,varargin)
	% make a column vector
	if (isvector(x) && isrow(x))
		x = cvec(x);
		ir = true;
	else
		ir = false;
	end
	
	method = 'linear';
	n1     = size(x,1);
	n2     = size(x,2);
%	np     = length(p);
%	q      = NaN(np,nx);

	% expand w to size of x
	if (isvector(w))
		w = repmat(cvec(w),1,n2);
	end

	% for negative weights invert x and use absolute weigth
	x = bsxfun(@times,sign(w),x);

	[xs sdx] = sort(x,1);

	% reorder weights in order of x
	for idx=1:n2
		w(:,idx) = w(sdx(:,idx),idx);
	end
	% normalise cdf
	% shift by 1/2 is necessary to obtain unbiased estimate
	ws       = bsxfun(@times,cumsum(abs(w))-0.5*abs(w),1./sum(abs(w)));
	% TODO duplicates can occur if weights are zero, which throws an error during interpolation

	if (length(x) == 1)
		% constant extrapolation
		q = xs*ones(size(p));
	else
		for jdx=1:size(ws,2)
			fdx = [true; diff(ws) > sqrt(eps)];
			q(:,jdx)        = interp1(ws(fdx,jdx),xs(fdx,jdx),p,method,NaN);
		end
	end

	% restore row vectors
	if (ir)
		q = rvec(q);
	end
end % wquantile


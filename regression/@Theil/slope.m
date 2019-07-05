% Di 12. Jan 17:07:54 CET 2016
% Karl Kastner, Berlin
%% fit the slope with the Theil-Sen method
%
% robust estimation of the slope
function [slope, qs, obj] = slope(obj,X,Y,w,p)
	if (nargin() < 3)
		w = [];
	end
	n     = length(X);
	flag  = triu(true(n),1);
	XX    = repmat(X,1,n);
	dx    = XX-XX';
	for idx=1:size(Y,2)
	YY    = repmat(Y(:,idx),1,n);
	dy    = YY-YY';
	if (~obj.repeated_medians)
	dy_dx = dy(flag(:))./dx(flag(:));

	if (~isempty(w))
		% average of both weights
		W = 0.5*(bsxfun(@plus,cvec(w),rvec(w)));
		W = W(flag)/sum(W(flag));
		if (nargout() < 2)
			slope = wmedian(W,dy_dx);
		else
			qs    = wquantile(W,dy_dx,p);
			slope = qs(2);
		end
	else
		if (nargout() < 2)
			slope(idx) = median(dy_dx);
		else
			qs(:,idx)    = quantile(dy_dx,p);
			slope(idx) = qs(2,idx);
		end
	end % isempty w
	else % repeated medians
		dy_dx(1:n+1:n^2) = dy_dx(end,1:end);
		slope(idx)       = median(median(dy_dx(1:end-1,:)));
	end % else of if ~repeated_medians
	end % for idx
end % slope


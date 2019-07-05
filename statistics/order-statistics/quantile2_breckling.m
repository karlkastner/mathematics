% 2016-03-01 15:26:30.465880523 +0100
% Karl Kastner, Berlin
%% qunatile regression
function Q = quantile2_breckling(X,p)
	n = size(X,2);
	me   = median(X,2);
	dX   = bsxfun(@minus,X,me);
	if (nargin() < 3 || isempty(u))
		dir = dX;
	end
	dir  = bsxfun(@times,dX,1./pnorm(X));
	for idx=1:n 
		q = [0; 0];
		q = lsqnonlin(@(q) objective(dX,p,q,dir(:,idx)).^2 ,q);
		Q(:,idx) = me + q;
	end
	
end

function f = objective(X,p,Q,dir)
	X_min_Q = bsxfun(@minus,X,Q);
	%cosa    = (Q'*X_min_Q)./(pnorm(Q).*pnorm(X_min_Q));
	cosa    = (dir'*X_min_Q)./(pnorm(X_min_Q));
%	cosa = abs(cosa);
	w   = 1 - (1-2*p)*cosa;
%	c   = 0.01;
	c = 0;
	fdx = pnorm(X_min_Q) < c;
	w(fdx)  = w(fdx)/c;
	w(~fdx) = w(~fdx)./pnorm(X_min_Q(:,~fdx));

%	sum(w < 0)
	wX_min_Q = bsxfun(@times,X_min_Q,w);
	f = sum(wX_min_Q, 2);
end

function p = pnorm(X)
	p = sqrt(sum(X.^2,1));
end

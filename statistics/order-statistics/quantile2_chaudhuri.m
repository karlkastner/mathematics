% 2016-03-01 10:09:24.480878731 +0100
% Karl Kastner, Berlin
%% quantile regression
function Q = quantile_chaudri(X,p,u)
	if (nargin() < 3 || isempty(u))
		u = bsxfun(@minus,X,median(X,2));
	end
	nu  = sqrt(sum(u.^2));
	fdx = find(nu > 0);
	if (~isempty(fdx))
		u(:,fdx) = bsxfun(@times,u(:,fdx),1./nu(fdx));
	end
	u = (1-2*p)*u;
	Q = zeros(size(u));
%	mu = mean(X,2);
	me = median(X')';
%	X_ = bsxfun(@minus,X,me);
	for idx=1:size(u,2)
		ret = check_degeneracy(X,u(:,idx));
		if (ret > 0)
			Q(:,idx) = X(:,ret);
		else
			%Q(:,idx) = lsqnonlin(@(q) objective(X,q,u(:,idx)),me);
			Q(:,idx) = optimise(X,u(:,idx),me);
		end
	end
end

function ret = check_degeneracy(X,u)
	n = size(X,2);
	ret = 0;
	for idx=1:n
		dX = bsxfun(@minus,X,X(:,idx));
		ndX = 1./pnorm(dX);
		f = sum(bsxfun(@times,dX,1./ndX),2) - (n-1)*u;
		if (pnorm(f) < 1 + pnorm(u))
			ret = idx;
			return;
		end
	end
end

function Q = optimise(X,u,Q)
	n = size(X,2);
	m = size(X,1);
	I = eye(m);
	for idx=1:10
		X_min_Q = bsxfun(@minus,X,Q);
		nX_min_Q = pnorm(X_min_Q);
		Phi = zeros(m);
		for jdx=1:n
			Phi = Phi + 1./nX_min_Q(jdx)*(I - 1./nX_min_Q(jdx).^2*X_min_Q(:,jdx)*X_min_Q(:,jdx)');
		end
		delta = sum( bsxfun(@times,X_min_Q,1./nX_min_Q),2) + n*u;
		
		Qold = Q;
		Q = Q + Phi \ delta;
		if (sum((Qold - Q).^2) < 1e-3)
			break
		end
	end
end
	
function f = objective(X,Q,u)
	X_min_Q = bsxfun(@minus,X,Q);
	nX_min_Q = pnorm(X_min_Q);
	n = size(X,2);
	f = sum( bsxfun(@times,X_min_Q,1./nX_min_Q),2) + n*u;
end

function f = pnorm(a)
	p = 2;
	f = sum(abs(a).^p,1).^(1/p);
end


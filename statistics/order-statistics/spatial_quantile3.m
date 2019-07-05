% Mo 18. Jan 09:02:10 CET 2016
% Karl Kastner, Berlin
%% spatial quantile
function q = spatial_quantile3(X,a,u)
	MaxIter = 20;
	abstol = 1e-7;
	% initial condition
	q = mean(X,2);
	if (0)
	q = lsqnonlin(@(q) func(X,u,q),q);
	else
	% iteration
	idx=0;
	while (true)
		[delta, Phi] = func(X,u,q);	
		q_ = q;
		q = q + Phi \ delta;
		d = norm(q-q_);
		if (d < abstol)
			break;
		end
		idx=idx+1;
		if (idx == MaxIter)
			m = size(X,1);
			q = NaN(m,1);
			break;
		end
	end
	idx
	end
end

function [delta Phi] = func(X,u,q)
	m = size(X,1);
	n = size(X,2);
	I = eye(m);
	p = 2;
	X_min_Q = bsxfun(@minus,X,q);	
	b = 1./lpnorm(X_min_Q,p);
	delta   = sum(bsxfun(@times,b,X_min_Q),2) + n*u;
	if (nargout > 1)
	if (0)
		Phi = zeros(m);
		for jdx=1:n
			Phi     = Phi + b(jdx)*(I - b(jdx).^2*X_min_Q(:,jdx)*X_min_Q(:,jdx)');
		end
	else
		A(:,1,:) = X_min_Q;
		B(1,:,:) = X_min_Q;
		b_(1,1,:) = b;
		Phi = bsxfun(@times,b_,bsxfun(@minus,I,bsxfun(@times,b_.*b_,bsxfun(@times,A,B))));
		Phi = squeeze(sum(Phi,3));
	end
	end
end


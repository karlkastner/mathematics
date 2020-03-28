% Fri Jun  8 15:46:25 MSK 2012
% Karl Kästner, Berlin

% V'VT is DT, where D is diag(D)

%recycling - block not possible
function [x V T r beta R R2] = minres_recycle(A, b, tol, V_old, T_old, maxiter)
	n = size(A,1);
	% preallocate memory
	V = zeros(n,maxiter);
	beta = zeros(1,maxiter);
	alpha = zeros(1,maxiter);
	%T = [];
	%V = [];
	% initial guess
	x = rand(n,1);
%	x = b;
	% non-trivial start residual
	r = b - A*x;
	% recycle old information
	if (~isempty(V_old))
		'Honk'
		%r = r - (r'*V)*V;
		x = x + V_old*(T_old\(V_old'*r));
		r = A*x;
	end
	beta(1) = sqrt(r'*r);
	V(:,1)  = r;
	r_ = zeros(size(r));
	r_(1) = beta(1);
	r0_ = r_;
	Q_=eye(length(b));
	r__ = r_;

%	beta(1) = sqrt(b'*b);
%	V(:,1)  = b;

	% start the new lanczos iteration
	% no restart like in gmres as matrix is symmetric
	idx=1;
	while (1)
		V(:,idx) = (1/beta(idx))*V(:,idx);
		if (idx > 1)
			V(:,idx+1) = A*V(:,idx) - beta(idx)*V(:,idx-1);
		else
			V(:,idx+1) = A*V(:,idx);
		end
		alpha(idx)  = V(:,idx+1)'*V(:,idx);
		V(:,idx+1)  = V(:,idx+1) - alpha(idx)*V(:,idx);
		beta(idx+1) = sqrt(V(:,idx+1)'*V(:,idx+1));
		% full reorthogonalisation
%		if (idx > 2)
%			% gram schmidt
%			%V(:,idx+1) = V(:,idx+1) - V(:,1:idx-2)*(V(:,idx+1)'*V(:,1:idx-2))';
%			% modified gram schmidt
%			for jdx=1:idx-2
%				V(:,idx+1) = V(:,idx+1) - (V(:,idx+1)'*V(:,jdx))*V(:,idx+1);
%			end
%			
%			V(:,idx+1) = beta(idx+1)/(V(:,idx+1)'*V(:,idx+1))*V(:,idx+1);
%		end
		% compute the true residual
		% construct the T matrix
%		end
		if (idx > 1)
			T = diag(sparse(beta(2:idx)),-1) + diag(sparse(alpha(1:idx))) + diag(sparse(beta(2:idx)),+1);
			% construct the final solution
			x = x + V(:,1:idx)*(T\(beta(1)*[1; zeros(idx-1,1)]));
			R(idx) = norm(A*x - b);
			R2(idx) = norm(A*x - (A\b));
			%R(idx) = norm(A*x - b)/norm(x);
		end
		D(idx) = alpha(idx);
		% perform a given's transform to make T tridiagonal
		if (idx > 1)
			% entries of the 2x2 given's block [c s; -s c]
			den   = sqrt(D(idx-1)^2  + beta(idx)^2);
			s     =  beta(idx)/den;
			c     =  alpha(idx-1)/den;
			% update the residual
			r_1 =  c*r_(idx-1);
			r_2   = -s*r_(idx-1);
			r_(idx-1) = r_1;
			r_(idx)   = r_2;
			d_1   = -(D(idx-1)^2      + beta(idx)^2)/den;
			d_2   = (D(idx-1)*D(idx) - beta(idx)^2)/den;
			e_2   = -(D(idx-1)*beta(idx) + beta(idx)*D(idx))/den;
			D(idx-1) = d_1;
			D(idx)   = d_2;
			E(idx)   = e_2;
			D_ = diag(sparse(D(1:idx))) + diag(sparse(E(2:idx)),+1);
			T = diag(sparse(beta(2:idx)),-1) + diag(sparse(alpha(1:idx))) + diag(sparse(beta(2:idx)),+1);
			full(D_)
			%full(T)
			T_ = T;
			T_(length(r__),length(r__)) = 0;
			[q_ r___] = qr(full(T_));
			Q =  eye(length(r));
			Q(idx-1,idx-1) = c; Q(idx,idx)   = c;
			Q(idx-1,idx) = s;   Q(idx,idx-1) = -s;
			Q_ = Q_*Q;
			%full(q_)
			%Q_
			%Q = q_; %*Q;
			%norm(A\b - x)/norm(x)
			%[(A\b - x) r__]
			r__ = q_*r__;
			[r_ r__ A*x-b]
			pause
		end
		%s(idx) = alpha(idx)/sqrt(h11^2 + h21^2);
		%c(idx) = alpha(idx)/sqrt(h11^2 + h21^2);

%		S = 

		if (beta(idx+1) < tol)
			break;
		end
		if (idx+1 > maxiter)
			'no convergence'
			break;
		end
		idx=idx+1;
	end
	% construct the T matrix
	T = diag(sparse(beta(2:idx)),-1) + diag(sparse(alpha(1:idx))) + diag(sparse(beta(2:idx)),+1);
	% construct the final solution
	x = x + V(:,1:idx)*(T\(beta(1)*[1; zeros(idx-1,1)]));
	r = V(:,idx+1);
	V = V(:,1:idx);
end % function minres


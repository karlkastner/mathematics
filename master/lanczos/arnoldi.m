% Fri Nov 19 14:36:40 MSK 2010
% Karl Kästner, KTH Stockholm

% incomplete hessian/qr-factorisation
% A*Q(:,1:end-1) - Q*H = 0

function [Q H] = arnoldi(A,b,n)
	l = size(A,1);
	Q = zeros(l,n);
	H = zeros(l,l); % l,n
	if (isempty(b))
		b = ones(l,1); % arbitrary ?
	end
	Q(:,1) = b / norm(b);
	r_old = 0;
	for idx=1:n
		% construct H = Q'*A*Q
		v = A*Q(:,idx);
		for jdx=1:idx
			H(jdx,idx) = Q(:,jdx)'*v;
			v = v - H(jdx,idx)*Q(:,jdx);
		end
		% v = A*Q(:,1:end-1) - Q*H (residual)
		% v_old'*v = 0
		nv = sqrt(v'*v);
		H(jdx+1,jdx) = nv;
		Q(:,jdx+1) = v / nv;
		% break down
		if (nv < 1E-12)
			H=H(1:idx+1,1:idx);
			Q=Q(:,1:idx+1);
			return
		end
	end
	H = H(1:idx+1,1:idx);
end % function arnoldi


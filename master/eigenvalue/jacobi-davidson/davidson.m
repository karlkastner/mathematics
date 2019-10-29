% Sun Jan  8 01:16:03 MSK 2012
% Karl Kästner, Berlin

%  n=10; A = n^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n), [v e] = eig(full(A)); [V E k] = davidson(A), diag(e)'
%	https://vpn.lan.kth.se
 
% Davidson Method
function [V_ E R] = davidson(A, B)
	abstol = sqrt(eps);;
	n = size(A,1);
	if (nargin < 2)
		B = speye(n);
	else
		G = zeros(n);
		gehep = 1;
	end
	I = speye(n);
	D = diag(sparse(diag(A)));
	%V = rand(n,1);
	% guess initial search space
	V = ones(n,1);
	V=V/norm(V);
	k = 1;
	H = zeros(n,n);
	while (1)
		% expand projected system
		AV(:,k) = A*V(:,k);
		H(1:k,k) = V(:,k)'*AV;
		H(k,1:k) = H(1:k,k)';
		if (ghep)
			BV(:,k) = B*V(:,k);
			G(1:k,k) = V(:,k)*BV;
			G(k,1:k) = G(1:k,k)
		end
		% find approximate eigenvalue
		if (~ghep)
			[y theta] = eigs(H(1:k,1:k), [], 1,'LA');
		else
			[y theta] = eigs(H(1:k,1:k), G(1:k,1:k), 1,'LA');
		end
		% approximated eigenvector
		x = V*y;
 		% residual (r = A*x - theta*B*x)
		if (~ghep)
			r = AV*y - theta*x;
		else
			r = AV*y - theta*BV*y;
		end
		% check for convergence
		R(k) = norm(r);
		if (R(k) < abstol)
			E(1) = theta;
			V_(:,1) = V(:,end);
			return;
		end
		% expand search space
		V(:,k+1) = gmres(D - theta*I, r);
		% orthogonalise by modified gram schmidt
		V(:,k+1) = mgs(V(:,1:k), V(:,k+1), k);
		k = k+1;
	end % while 1
end % davidson


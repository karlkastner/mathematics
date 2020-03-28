% Mon Oct 24 18:50:54 MSK 2011
% Karl Kästner, Berlin

function [X L K R] = jacobi_davidson_simple(A)
	MaxIter = 100;
	tol = 1e-7;
	n = size(A,1);
	t = ones(n,1);
	R = [];
	X = [];
	L = [];
	M =  [];
	V =  [];
	AV = [];
	m=1;
	while(1) % m
		% orthogonalise vector to current subspace
		t = mgs(V,t,m-1);
		V(:,m) = t;
		% complete stored matrix-subspace product
		AV(:,m) = A*V(:,m);
		% latest update of projected system M = V'*A*V
		M(1:m,m) = V(:,m)'*AV;
		M(m,1:m) = M(1:m,m)';
		% find eigenpair of projected system
		[s theta] = eigs(M(1:m,1:m),[],1,'SA');

theta
pause
		% ritz-vector (approximated eigenvector)
		u  = V*s;
		% residual, AV*s = A*u
		r  = AV*s - theta*u;
		% check for convergence
		R(m) = norm(r);
		if (R(m) < tol)
			L = theta;
			X = u;
			return;
		end
		% solve: ((I - u*u')*(A - theta*I)*(I - u*u'))t = -r
		% rank deficient matrix, iterative solution obligatory
		t = gmres(@(x) afun_jdm(A,u,theta,x), -r);
	
		if (MaxIter == m)
			'error no convergence'
			return
		end
		m=m+1;
	end % while (1) % m 
end % jacobi_davidson


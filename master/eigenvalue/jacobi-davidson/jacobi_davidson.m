% Mon Oct 24 18:50:54 MSK 2011
% Karl Kästner, Berlin

function [X L K] = jacobi_davidson(A)
	MaxIter = 100;
	tol = 1e-7;
	X = [];
	L = [];
	
	nn = size(A,1);
	K = zeros(1,nn);
	ndx=1;
	while (1)
	%for ndx=1:nn
	n = nn-ndx+1;
	t = ones(nn,1);
	M = [];
	V = [];
	VA = [];
	mdx=0;
	while(1)
		% todo, repeat, if necessary
		for idx=1:mdx-1
			t = t - (V(:,idx)'*t)*V(:,idx);
		end % for idx
		mdx=mdx+1;
		V(:,mdx) = t/norm(t);
		VA(:,mdx) = A*V(:,mdx);
		% latest update for M = V'*A*V
		for idx=1:mdx
			M(idx,mdx) = V(:,idx)'*VA(:,mdx);
			M(mdx,idx) = M(idx,mdx);
		end % for idx
		% todo, use symmetrie in M
		 %[s theta] = eigs(M,1,'LM')
		 %[s theta] = eigs(M,[],1,'LA');
		 % u=V*s;
		 %Au = VA*s; % equals uA = A*u
		 %r = Au - theta*u;
		[S Theta] = eig(M) % sorted pairs
		E = diag(Theta);
		[E pos] = sort(E);
		S = S(:,pos);
		Theta = diag(E);

		S
		Theta

		u = V*S(:,1);
		Au = VA*S(:,1);
		r = Au - Theta(1,1)*u;
		% check for convergence
		 %if (normr < tol) 
		while (norm(r) < tol)
			L(ndx,ndx) = Theta(1);
			X = [X u];
			ndx = ndx+1;
			if (ndx > nn) return; end;
			mdx=mdx-1;
			M = 0;
			for idx=1:mdx
				V(:,idx) = V*S(:,idx+1);
				VA(:,idx+1) = VA*S(:,idx+1);
				M(idx,idx) = Theta(idx+1,idx+1);
				Theta(idx,idx) = Theta(idx+1,idx+1);
				S(:,idx) = 0; S(idx,idx) = 1;
			end
			u = V(:,1);
			r = VA(:,1) - Theta(1,1)*u;
		end
		mmax=10;
		mmin=2;
		if (mdx > mmax)
			M = 0;
			for idx=2:mmin
				V(:,idx) = V*S(:,idx);
				VA(:,idx) = VA*S(:,idx);
				M(idx,idx) = Theta(idx,idx);
			end
			V(:,1) = u;
			VA(:,1) = Au;
			M(1,1) = Theta(1,1);
			mdx=mmin;
		end

		I = eye(nn);
		% todo, approximate solver
		Q = [X u];
		%C = ((I - u*u')*(A - theta*I)*(I - u*u'));
		Theta = Theta(1,1)
		C = ((I - Q*Q')*(A - Theta*I)*(I - Q*Q'));
		%t = C \ -r; % rank deficient ?
		t = gmres(C, -r);
	
		% test orthogonality
		norm(t'*Q); %norm(t'*u)
	
		K(ndx)=K(ndx)+1;
		if (MaxIter == mdx)
			'error no convergence'
			return
		end
	end % while (1)
	
	% deflate, attention, loss of sparsity !!!
	% orthonormal basis
	%B = [X(:,1:ndx) [zeros(1,n-1); eye(n-1)]];
%	B = [u [zeros(1,n-1); eye(n-1)]];
	%for idx=ndx+1:n
%	for idx=2:n
%	        B(:,idx) = B(:,idx) - B(:,1)'*B(:,idx)*B(:,1);
		%for jdx=1:ndx
	        %	B(:,idx) = B(:,idx) - B(:,jdx)'*B(:,idx)*B(:,jdx);
		%end % for jdx
%	end % for idx
%	sort(eig(A))
	%A = B'*A*B
	%A = inv(B)*A*B
%	A = A(2:end,2:end);
	%A = A(ndx+1:n,ndx+1:n);
%	sort(eig(A))
	end % for ndx
end % jacobi_davidson


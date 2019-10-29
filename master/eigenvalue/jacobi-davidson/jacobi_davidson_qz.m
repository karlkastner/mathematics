% Thu Jan 26 14:14:32 MSK 2012
% Karl Kästner, Berlin

function [Q Z Lambda] = jacobi_davidson_qz(A,B,k_max,SIGMA)
	tol = sqrt(eps);
	n = size(A,1);
	maxiter = n;
	solver_maxiter = n;
	m_max = 1000;%k_max + 10; %100;
	m_min = 5;

	t = [];
	Lambda = [];
	V  = [];
	AV = [];
	BV = [];
	Q  = [];
	Z  = [];
	if (isempty(t))
		t = ones(n,1) + 0.1*rand(n,1);
		t = t/norm(t);
	end

	m = 1;
	m_ = 1;
	k = 1;
	while (k <= k_max)
		[V(:,m) BV(:,m)] = mgs(V,BV,B,t,m-1);
		AV(:,m)  = A*V(:,m);
		M(1:m,m) = V(:,m)'*AV;
		M(m,1:m) = M(:,m)';

		[S Theta] = eig(M);
		Theta = diag(Theta);
		switch (SIGMA)
			case {'LM'}
				[void p] = sort(abs(Theta),1,'descend');
			case {'SM'}
				[void p] = sort(abs(Theta),1,'ascend');
			case {'LA'}
				[void p] = sort(Theta,1,'descend');
			case {'SA'}
				[void p] = sort(Theta,1,'ascend');
		end
		Theta = Theta(p);
		S = S(:,p);
		u = V*S(:,1);
		Q(:,k) = u;
		Z(:,k) = BV*S(:,1);
		r = AV*S(:,1) - Theta(1)*Z(:,k);
		while (norm(r)/norm(u) < tol)
			Lambda(k,1) = Theta(1,1);
			if (k == k_max)
				if (nargout < 2)
					Q = Lambda;
				else
					Lambda = diag(Lambda);
				end
				return;
			end
			k = k+1;

			m = m-1;
		
			V  = V*S(:,2:m+1);
			AV = AV*S(:,2:m+1);
			BV = BV*S(:,2:m+1);
			M  = diag(Theta(2:m+1,1));
			S  = eye(m);
			Theta = Theta(2:m+1,1);
		
			if (0 == size(V,2))
				'start vector was an eigenvector, unable to proceed'
				return
			end

			u = V(:,1);
			Q(:,k) = u;
			Z(:,k) = BV(:,1);
			r = AV(:,1) - Theta(1)*Z(:,k);
		end  % while norm(r)/norm(u) < tol
		% restart - delfate the m-matrix
		if (m >= m_max)
			V_old = V;
			AV_old = AV;
			BV_old = BV;
			for idx=2:m_min
				V(:,idx)  = V_old*S(:,idx);
				AV(:,idx) = AV_old*S(:,idx);
				BV(:,idx) = BV_old*S(:,idx);
				M(idx,idx) = Theta(idx,1);
			end
			V(:,1) = u;
			AV(:,1) = AV_old*S(:,1); 
			BV(:,1) = BV_old*S(:,1); 
			M(1,1) = Theta(1);
			m = m_min;
		end % restart

		v=ver;
		if (1 == strcmp(v(1).Name,'Octave')) % inverted c-logic
			% octave
			[t res] = pgmres(@(x) afun_jdm(A, X, Theta(1,1), x), -r, zeros(size(t)), tol, 10, 10, speye(n));
			if (res(end) > tol)
				flag = 1
			end
		else
			%matlab
		        %[t flag] = gmres(@(x) afun_jdm(A, B, Q, Z, Theta(1,1), x), -r, [], tol, maxiter, [], [], zeros(size(A,1),1));
		        [t flag] = minres(@(x) afun_jdm(A, B, Q, Z, Theta(1), x), -r, tol, solver_maxiter, [], [], zeros(size(A,1),1));
		end
		if (0 ~= flag)
			'error, solver minres stagnated'
			if (1 == nargout())
				Q = [Lambda Theta(1)];
			end
			return;
		end
		m  = m+1;
		m_ = m_+1;
		if (m_ > maxiter)
			'error: maximum number of iterations reached'
			if (1 == nargout())
				Q = [Lambda Theta(1,1)];
			else
				Lambda = diag(Lambda);
			end
			return;
		end
	end % while k <= k_max
end % function jacobi_davidson_qz


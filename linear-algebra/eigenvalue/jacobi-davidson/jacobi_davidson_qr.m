% Sun Jan 22 22:28:34 MSK 2012
% Karl Kästner, Berlin
%
function [X Lambda] = jacobi_davidson_qr(A, k_max, SIGMA, V0, tau)
	tol   = sqrt(eps);
	n = size(A,1);
	maxit = n;
	solver_maxit = n;
	solver_tol = 1e-1;
	m_min = k_max+5;
	m_max = m_min+5;
	version=ver();

	t = [];
	Lambda = []; 
	X = [];
	if (nargin()<4 || isempty(V0))
	        % build an initial search space (Arnoldi)
	        % this is not in sleipen 1999
		V0 = ones(n,1) + 0.1*rand(n,1);
	        m  = min(n, m_min); % max ???
	        [V H] = arnoldi(A, m-1, V0);
		m = min(m,size(V,2));
	else
		V = orth(V0);
		m = size(V,2)
	end
	AV = full(A*V);
	M = V'*AV;

	m_ = 1;
	k  = 1;

	while (k <= k_max)
		[S Theta] = eig(M);
		Theta = diag(Theta);
		% eigenvalue selection (added)
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
		X(:,k) = u;
		r = AV*S(:,1) - Theta(1)*u;
		while (norm(r)/norm(u) < tol) % added / norm(u)
			Lambda(k,1) = Theta(1,1);
			if (k == k_max)
				if (nargout < 2)
					X = Lambda;
				else
					Lambda = diag(Lambda);
				end
				return;
			end
			k = k+1;
			% shrink search space by converged eigenpair
			m = m-1;

			% the next two steps is incorrectly stated (Sleij 1999)
			V  = V*S(:,2:m+1);
			AV = AV*S(:,2:m+1);

			M  = diag(Theta(2:m+1,1));
			S  = eye(m);
			Theta = Theta(2:m+1,1);

			% the next step is missing in (Sleij 1999)
			if (0 == size(V,2))
				% search space empty, get a new start vector
				% TODO, get a new start space of rank j_min
				V  = rand(n,1); V = (1/sqrt(V'*V))*V;
				AV = A*V;
				M  = V'*AV;
				Theta = M;
			end
			u = V(:,1);
			X(:,k) = u;
			if (k == k_max) r = 0; end
			r = AV(:,1) - Theta(1)*u;
		end % while ||r|| < tol
		% deflation
		if (m >= m_max)
			% in case of convergence there is no deflation and
			% u is therefore always V(:,1) at deflation
			% V  := [u   V*S(:,2:m_min)];
			% AV := [A*u AV*S(:,2:m_min)];
			m = m_min;
			V  = V*S(:,1:m);
			AV = AV*S(:,1:m);
			M = diag(Theta(1:m));
		end % deflation

		if (1 == strcmp(version(1).Name,'Octave')) % inverted c-logic
			% octave
			[t res] = pgmres(@(x) afun_jdm(A, X, Theta(1,1), x), -r, zeros(size(t)), tol, 10, 10, speye(n));
			if (res(end) > tol)
				flag = 1
			end
		else
			% matlab
			% the important projection of the residual is missing in Sleij 1999
			r_ = r - X*(X'*r);
			% todo, keep shift tau for the first iteration steps konstant
			%if (norm(r) < 1e-3)
				tau = Theta(1);
			%end
			x0 = zeros(size(u));
			%x0 = u; x0 = x0 - X*(X'*x0); % makes no difference
			[t flag] = minres(@(x) afun_jdm(A, [], X, X, tau, x), -r_, solver_tol, solver_maxit, [], [], x0);
			%t  = t - X*(X'*t); % makes also no difference
			%[t flag] = gmres(@(x) afun_jdm(A, [], X, X, Theta(1), x), -r, [], solver_tol, solver_maxit, [], [], zeros(size(r)));
		end
		if (0 ~= flag)
			'error, solver stagnated'
			[norm(r) norm(u) norm(r)/norm(u) norm(A*V(:,1) - AV(:,1)) norm(V'*V-eye(m)) norm(X'*X - eye(k))]
			if (nargout() < 2)
				X = [Lambda; Theta(1,1)];
			else
				Lambda = diag(Lambda);
			end
			return
		end
		m = m+1;
		m_ = m_+1;
		if (m_ > maxit)
			'error: maximum number of iterations reached'
			if (nargout() < 2)
				X = [Lambda; Theta(1,1)];
			else
				Lambda = diag(Lambda);
			end
			return
		end
		% expand the search space
		% shifted from beginning to end of loop
		V(:,m) = mgs(V,[],[],t,m-1); % mgs returns normed result
		AV(:,m) = A*V(:,m);
		M(1:m,m) = V(:,m)'*AV;
		M(m,1:m) = M(:,m)';	% added
	end % while k < k_max
end % function jacobi_davidson_qr


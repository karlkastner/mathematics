% Wed Jun 13 01:20:05 MSK 2012
% Karl Kästner, Berlin

% A=rand(10); [e v] = inverse_iteration(A,[],[],[],[]); [v_ e_] = eigs(sparse(A),[],1); e-e, norm(v.^2-v_.^2)

function [lambda x R] = inverse_iteration(A, x, tol, maxiter, flag)
	if (isempty(tol))
		tol = 1e-12;
	end
	if (isempty(maxiter))
		maxiter = 100;
	end
	if (isempty(x))
		x = rand(size(A,1),1);
	end
	if (isempty(flag))
		flag = 0;
	end
	I = speye(size(A));
	x = 1/sqrt(x'*x)*x;
	lambda = 0;
	idx=0;
	% precompute LU-factorisation of A
	if (2 == flag)
                [L_,U_,pp_,qq_,D_] = lu(A);
	end
	while (1)
		x_old = x;
		switch(flag)
			case {0} % shift and exact
				sigma = lambda;
				x = (A - sigma*I) \ x;
			case {1} % shift and not exact
				x0=x;
				sigma = lambda;
				[x flag_ relres iter1] = minres(A-sigma*I,x,tol,maxiter,[],[],x0);
			case {2} % no shift and exact
				sigma=0;
                		x = qq_*(U_ \ (L_ \ (pp_*(D_\x))));
			case {3} % no shift and no exact
				sigma=0;
				x0 = x;
				[x flag_ relres iter1] = minres(A,x,tol,maxiter,[],[],x0);
		end
		lambda_old = lambda;
		lambda = 1/(x'*x_old) + sigma;
		x = 1/sqrt(x'*x)*x;
		idx = idx+1;
		R(idx) = abs(lambda - lambda_old);
		if (R(idx) < tol || idx > maxiter)
			%lambda = x'*(A*x); % get correct sign of lambda
			return;
		end
	end
end


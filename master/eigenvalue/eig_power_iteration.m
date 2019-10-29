% Wed Jun 13 00:45:18 MSK 2012
% Karl Kästner, Berlin

% A=rand(10); [e v] = power_iteration(A,[],[],[]); [v_ e_] = eigs(sparse(A),[],1); e-e, norm(v.^2-v_.^2)

function [lambda x L] = power_iteration(A, x, tol, maxiter)
	if (nargin() < 3 || isempty(tol))
		tol = sqrt(eps);
	end
	if (nargin() < 4 || isempty(maxiter))
		maxiter = 100*size(A,1);
	end
	if (nargin() < 2 || isempty(x))
		x = rand(size(A,1),1);
	end
	x = 1/sqrt(x'*x)*x;
	lambda = 1;
	idx=0;
	while (1)
		%x_old = x;
		x = A*x;
		lambda_old = lambda;
		lambda = sqrt(x'*x);
		x = (1/lambda)*x;
		%lambda = x'*x_old;
		%x = 1/sqrt(x'*x)*x;
		idx = idx+1;
		R(idx) = abs(lambda - lambda_old);
		R(idx)
		if ( R(idx) < tol || idx > maxiter)
			% get signed eigenvalue
			x_old = x;
			x = A*x;
			lambda = x'*x_old;
			x = (1/lambda)*x;
			return;
		end
	end % while 1
end % power_iteration


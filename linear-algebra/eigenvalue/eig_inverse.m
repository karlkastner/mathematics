% Sun Dec  2 15:02:12 MSK 2012
% Karl Kästner

% inverse iteration
function [x lambda R L] eig_inverse(A, x, lambda, tol)
	if (isempty(tol))
		tol = 1e-12;
	end
	if (isempty(x0))
		x = rand(size(A,1),1);
	end
	if (isempty(lambda0))
		lamda = 0;
	end
	I = speye(size(A));
	
	while(1)
		% x = (A - lambda*I) \ x;
		% todo, make A - lambda*I a function)
		% todo, tol
		x = minres(A - lambda*I, x);
		x = x / norm(x);
		Ax = A*x;
		lambda = x'*Ax;
		resnorm = norm(lambda*x - Ax);
		L(end+1) = lambda;
		R(end+1) = resnorm;
		if (resnorm < tol)
			break;
		end
	end
end % function eig_inverse


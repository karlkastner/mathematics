% Thu Jan 26 20:50:42 MSK 2012
% Karl Kästner, Berlin

function [Q T BQ] = lanczos_ghep(A,B,k)
	n = length(A);
	Q(:,1) = rand(n,1);
	BQ(:,1) = B*Q(:,1);
	beta(1) = sqrt(Q(:,1)'*BQ(:,1));
	for idx=1:k
		BQ(:,idx)   = (1/beta(idx))*BQ(:,idx);
		Q(:,idx)    = (1/beta(idx))*Q(:,idx);
		BQ(:,idx+1) = A*Q(:,idx);
		if (idx > 1)
			BQ(:,idx+1) = BQ(:,idx+1) - beta(idx)*BQ(:,idx-1);
		end
		alpha(idx)  = Q(:,idx)'*BQ(:,idx+1);
		BQ(:,idx+1) = BQ(:,idx+1) - alpha(idx)*BQ(:,idx);
		%Q(:,idx+1)  = B \ BQ(:,idx+1); % todo, use cg / minres
		[Q(:,idx+1) flag] = minres(B, BQ(:,idx+1));
		if ( 0 ~= flag)
			'minres did not converge'
			return
		end
		beta(idx+1) = sqrt(abs(Q(:,idx+1)'*BQ(:,idx+1)));
	end % for idx
	T =	  diag(sparse(beta(2:end-1)),-1) ...
		+ diag(sparse(alpha)) ...
		+ diag(sparse(beta(2:end-1)),+1);
	Q  = Q(:,1:k);
	BQ = BQ(:,1:k);
end % lanczos_ghep


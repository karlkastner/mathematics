% Mon Jan  2 18:57:51 MSK 2012
% Karl Kästner, Berlin

% Test : 
%> A =rand(10,5); B=rand(10); B=B*B'; L=chol(B); V = inv(L)*orth(A); x = rand(10,1); I = eye(10);
%>> y = (I - V*V'*B)*x, y=y/sqrt(y'*B*y), y'*B*[V y]             
% modified gram schmidt algorithm
function [x Bx] = mgs(Q, BQ, B, x, n)
	% todo, repeat, if necessary
	if (isempty(B))
		n0 = norm(x);
		for idx=1:n
			x = x - (Q(:,idx)'*x)*Q(:,idx);
		end % for idx
		n1 = norm(x);
			%if (n1 < n0*0.5)
			while (n1 < n0*0.5)
				n0 = n1;
				for idx=1:n
					x = x - (Q(:,idx)'*x)*Q(:,idx);
				end % for idx
				n1 = norm(x);
			end
		x = x/n1;
	else % generalised inner product
		for idx=1:n
			x = x - (BQ(:,idx)'*x)*Q(:,idx);
		end
		Bx = B*x;
		inBx = 1/sqrt(x'*Bx);
		x  = inBx*x;
		Bx = inBx*Bx;
	end
end % mgs


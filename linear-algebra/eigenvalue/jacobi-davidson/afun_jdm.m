% indirectly computes the matrix vector product
function x = afun_jdm(A, B, u, v, theta, x) %, L)
	%C = ((I - u*u')*(A - theta*I)*(I - u*u'));
	%x = C*x;
%	x = x - (u'*x)*u;
%	x = A*x - theta*x;
%	x = x - (u'*x)*u;
	x = x - u*(v'*x);
	if (~isempty(B))
		x = A*x - theta*B*x;
	else % B = I
		% precondition
%		if (nargin() > 6 && ~isempty(L))
%			x = L'\(L\(A*x)) - theta*x;
%		else
			% no preconditioning, L = I
			x = A*x - theta*x;
%		end
	end
	x = x - v*(u'*x);
end


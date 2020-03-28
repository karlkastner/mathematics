% 2012-04-19 13:02:01 UTC
% Karl Kästner, Berlin

function derive_fdm_error()
	n = 2;
	m =4;	
	syms h
	syms d1f d2f d3f d4f
	dF = [d1f d2f d4f d4f].'

%	A = (sym([-n:-1 1:n].'*ones(1,m)).^(ones(2*n,1)*(1:m)))./(ones(2*n,1)*factorial(1:m)) * diag(h.^(1:m))
%	B = [0 -1 0; 1 -1 0; 0 -1 1; 0 -1 0]
	A = [-h  h^2/2;
              h  h^2/2]
	B = [ 1 -1 0    h^3/6 -h^4/24;
	      0 -1 1   -h^3/6 -h^4/24]
	
	K = inv(A)*B
	%K = A \ B %inv(A'*A)*(A'*B)
	%A.'*A	
%	B'*B
%	(B'*B) \ (B'*A)
%	C = (A'*A) \ (A'*B)
%	B'*A
%	(B'*A) \ (B'*B)
%	A*C
%	A*C - B
end % derive_fdm_error()


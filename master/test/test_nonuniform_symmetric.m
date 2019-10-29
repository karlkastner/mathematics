function test_nonuniform_symmetric()

n=4;

% 1D symmetric matrix A and scaling factor B
I=eye(4);
A = full(spdiags(ones(n,1)*[1 -2 1],-1:1,n,n));
B=0.5*eye(n);

x_ = (0:n+1)';
[D L] = laplacian_non_uniform(x_);
'1D test 0 is pass'
norm(L - L') + norm(D - B) + norm(A - L)


% 2D asymmetric and symmetric system
M = kron(B*A,I) + kron(I,B*A);
AA = (kron(I,inv(B))*kron(A,I) + kron(inv(B),I)*kron(I,A));
BB = kron(I,B)*kron(B,I);

[DD LL] = laplacian_combine(D,L,D,L);
'2D test symmetrie and correctness '
norm( LL - LL') + norm((AA)' - (AA)) + norm(M - (DD*LL))

%n=5;
%x = (0:n+1)'
%L = laplacian_non_uniform(x);
%full(L)

end % test_nonuniform_symmetric

function [DD LL] = laplacian_combine(D1,L1,D2,L2)
	I1 = speye(size(D1));
	I2 = speye(size(D2));
	DD = kron(D1,I2)*kron(I1,D2);
	LL = kron(I1,inv(D2))*kron(L1,I2) + kron(inv(D1),I2)*kron(I1,L2);
end % function laplacian_combine



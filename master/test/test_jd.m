% 1D
n=250;
A = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); 
tic
[V E No Ni] = jacobi_davidson_qr(A,10);
E
[No Ni]
toc
try
jdqr(A,10)
catch
	'Error'
end

% 2D
n=10;
A1 = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
I1 = speye(n);
n=25;
A2 = spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
I2 = speye(n);

A = kron(I1,A2) + kron(A1,I2);

tic
%jacobi_davidson_qr(A,1);
[V E No Ni] = jacobi_davidson_qr(A,10);
E
[No Ni]
toc
try
jdqr(A,10)
catch
	'Error'
end

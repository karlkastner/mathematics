% Mon 13 Jun 13:01:08 CEST 2022

n = [5,5];
L = 2*[1,1]
%n(1),n(2)];
x = randn(prod(n),1);

[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,'dirichlet');
D2 = D2x+D2y;

p = 1;
y = laplacian_power(L,n,p,x);
y(:,2) = D2*x;
y(:,1)./y(:,2)

if (0)
% inverse
p = -1;
y = laplacian_power(L,n,p,x);
y(:,2) = D2 \ x;
y(:,1)./y(:,2)

% sqrt
p = -1/2;
y = laplacian_power(L,n,p,x);
y(:,2) = sqrtm(full(D2)) \ x;
y(:,1)./y(:,2)
end

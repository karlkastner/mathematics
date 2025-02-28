
%	n=10; 
%	D1 = downsampling_matrix(n);
%	% x = (0:n-1)'; x=x.*x'; A*x*A'

n = 10;
%A = upsampling_matrix(n);

A = downsampling_matrix(n)
B = upsampling_matrix(n/2)

B*A
A*B

%x=(0.5:2:2*n-0.5)';
x = (0:n-1)'
A*x
x(end) = -1;
A*x

%x=(0.5:2:n-0.5)';
x = (0:2:n-2)';
%x(end)=-1.5; A*x,
%xr(1) = 2*n+0.5; A*xr
B*x
x(1) = 10;
B*x

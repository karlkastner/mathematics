% 2015-08-09 13:03:54.674520831 +0200

rho = 0.95;
n = 600;
a = acfar1(rho,n);
s = 1;
w = 1./(n - (0:n-1)');

a_hat = a + s*w.*randn(n,1);

E = [];
m = 1000;
R = linspace(-2,2,m)';
for idx=1:m
	a = acfar1(R(idx),n);
	E(idx,1) = norm(a - a_hat);
end

semilogy(R,E)


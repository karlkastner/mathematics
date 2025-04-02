L = 1e3;
n = 1e6;
x = innerspace(0,L,n);
dx = x(2)-x(1);

a = [3,4,5];

mu = [];
for idx=1:length(a)
	f = powerpdf(x,a(idx));
	mu(idx,1) = powerpdf_mean(a(idx));
	mu(idx,2) = sum(f.*x)*dx;
end
mu

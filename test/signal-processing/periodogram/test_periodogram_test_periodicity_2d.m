% Tue 20 Dec 14:34:34 CET 2022
pp = 0:0.01:0.2;
pp = [0.05:0.05:0.95];

n = 100;
L = 1;
x = repmat(linspace(0,L,n)',1,n);
y = x';
fc = 20;

p = [];
nf = 5;
bmsk = false(n,n);
m = 10;
m1 = n; %/10;
m2 = 4;
bmsk(1:m1,1:m2) = true;
d = [];
mu = [];
for idx=1:length(pp)
	idx
	z = pp(idx)*sin(2*pi*x*fc) + (1-pp(idx))*randn(n,n);

	ns = 1000;
	[p(idx,1),stat] = periodogram_test_periodicity_2d(z, [L,L], nf, bmsk,[],[],ns);
%	d(idx,1) = stat.f.d1;
%	d(idx,2) = stat.f.d2;
%	mu(idx,1) = stat.mu;
%	mu(idx,2) = stat.sd;
	
	p1(idx,1) = stat.p1;
	[p(idx,2),stat] = periodogram_test_periodicity_2d(z, [L,L], nf, []);
	p1(idx,2) = stat.p1;
	%, frmin, frmax)
	p
	p1
%	-> use fisher_moment2param to approximate distribution
end
%d
%periodogram_test_periodicity_2d.m

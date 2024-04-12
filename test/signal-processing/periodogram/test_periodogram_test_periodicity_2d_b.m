% 2023-07-07 08:50:08.203620297 +0200
bmsk = [];
fmsk = [];
nf = 3.5;
n = 50;
L = 1;
[fx,fy,frr] = fourier_axis_2d([n,n],[n,n]);
fmsk = true(size(frr)) & (fx>=0);
%bmsk = true(n);
p = [];
ns = 100;
m = 1e3;
for idx=1:m
	b = rand(n);
	[isp,stat] = periodogram_test_periodicity_2d(b, [L,L], nf, bmsk, fmsk);
	p(idx,2) = stat.pn;
	%p(idx,1) = stat.pn;
	[pn,stat,out] = periodogram_test_periodicity_2d_old(b, nf, bmsk, fmsk, ns);
	p(idx,1) = pn;
end
mean(p<=0.05)

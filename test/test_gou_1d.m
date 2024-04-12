% Tue 28 Nov 15:45:52 CET 2023
geo = false;
dynamic_order = true;
mm = 8*128;
L = 1024/4*mm;
n = 1024*mm;
pmu = 0;
psd = 1;
theta = 0.35;
%theta = 0.01;
order0 = 40;
order = order0+1;
N = [];
mu = [];
sd = [];
r =[];
idx=0;
n0 = n;
while (L/n<32)
	order
	idx=idx+1;
	N(idx) = n;
	if (geo)
	[y,cov_] = geometric_ar1_1d_grid_cell_averaged_generate(pmu,psd,theta,L,n,order);
	else
	[y,cov_] = ar1_1d_grid_cell_averaged_generate(pmu,psd,theta,L,n,order);
	end
	if (1 == idx)
		x = cvec(y);
	end
	k = n/2;
	mu(1,idx) = mean(x(1:k));
	sd(1,idx) = std(x(1:k));
	mu(2,idx) = mean(y(1:k));
	sd(2,idx) = std(y(1:k));
	r(1,idx)  = corr(x(1:k-1),x(2:k));
	r(2,idx)  = corr(y(1:k-1),y(2:k));

	A = downsampling_matrix(n,'pairwise');
	x = A'*x;
	n = n/2;
if (dynamic_order)
	order = order0*n0/n+1;
	%order = 2*order;
end
end
L./N
mu
sd
r
sd(:,2:end)./sd(:,1:end-1)


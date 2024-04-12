% 2023-07-21 21:45:00.097435071 +0200

lmu_a = 2;
lsd_a = 3;
lmu_b = lmu_a;
lsd_b = lsd_a;
lr = 0.7;

corr_ = logn_corr(lr,lmu_a,lmu_b,lsd_a,lsd_b)

[mu,sd] = logn_param2moment(lmu_a,lsd_a)

[lmu,lsd,lr] = logn_moment2par_correlated(mu,sd,corr_)



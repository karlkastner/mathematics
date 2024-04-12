mu = 2;
sd = 3;
[lmu, lsd] = logn_moment2par(mu,sd);

x=lognrnd(lmu,lsd,1e6,1);
mean(x)
std(x)


mu = 2;
sd = 3;
x = linspace(-5,-4);

c = normcdf(x,mu,sd);
[mu_,sd_] = normcdf_fit(x,c)


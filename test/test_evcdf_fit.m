x = linspace(-3,10)';
mu = 2;
b  = 3;

c = evcdf(x,mu,b);
c(:,2) = 1-exp(-exp((x-mu)/b));

clf
subplot(2,2,1)
plot(x,c)

[mu_,b_] = gumbel_cdf_fit(x,c(:,1))

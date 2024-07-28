

a = 1/2;
b = 1/3;
L = 20;

mu = lognpdf_mean(a,b)
sd = lognpdf_std(a,b)

mu_ = quad(@(x) x.*lognpdf(x,a,b),0,L)
s2 = quad(@(x) (x-mu).^2.*lognpdf(x,a,b),0,L)
sd_ = sqrt(s2)

mu - mu_
sd - sd_

figure(1)
n = 1e4;
x = innerspace(0,L,n);
plot(x,lognpdf(x,a,b))

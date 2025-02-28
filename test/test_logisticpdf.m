x = linspace(-3,3)';
mu = 2;
a  = 3;

p  = logistic_pdf(x,mu,a);
c  = logistic_cdf(x,mu,a);
x_ = logistic_inv(c,mu,a);

[mu_,a_] = logistic_cdf_fit(x,c)

figure(1);
clf
subplot(2,2,1)
plot(x,[p,cdiff(c)./cdiff(x)]); 

subplot(2,2,2)
plot(x,c)

subplot(2,2,3)
plot(x,x_)




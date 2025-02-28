x = linspace(-3,3)';
mu = 2;
a  = 3;
b  = 5;

p = generalized_logistic_pdf(x,mu,a,b);
c  = generalized_logistic_cdf(x,mu,a,b);
x_ = generalized_logistic_inv(c,mu,a,b);

figure(1);
clf
subplot(2,2,1)
plot(x,[p,cdiff(c)./cdiff(x)]); 

subplot(2,2,2)
plot(x,c)

subplot(2,2,3)
plot(x,x_)




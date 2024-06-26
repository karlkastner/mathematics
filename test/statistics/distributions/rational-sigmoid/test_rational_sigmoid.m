% 2024-05-17 09:03:25.401783503 +0200
% Karl Kastner, Berlin

a = [0.5,1,2];
k = [0.5,1,2];
mu = 0;

x = linspace(-5,5)';
c = rational_sigmoid_cdf(x,mu,a,2);
p = rational_sigmoid_pdf(x,mu,a,2);

clf
subplot(2,2,1)
plot(x,c);

subplot(2,2,3)
plot(x,p);
hold on
set(gca,'colororderindex',1)
plot(mid(x),diff(c)./diff(x),'--');

c = rational_sigmoid_cdf(x,mu,1,k);
p = rational_sigmoid_pdf(x,mu,1,k);
subplot(2,2,2);
plot(x,c)

subplot(2,2,4);
plot(x,p);
hold on
set(gca,'colororderindex',1)
plot(mid(x),diff(c)./diff(x),'--');

par = [1,2,3];
rng(0)
n = 100;
x = -10+20*rand(n,1);
c = rational_sigmoid_cdf(x,par(1),par(2),par(3))+1e-3*randn(n,1);

par = fit_rational_sigmoid_cdf(x,c);

x=2.2
c = rational_sigmoid_cdf(x,1,2,3);
rational_sigmoid_inv(c,1,2,3)


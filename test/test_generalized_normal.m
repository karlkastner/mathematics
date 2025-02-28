a = [0.75,1,1.5];
b = [0.75,1,1.5];
mu = 0;

x = linspace(-5,5)';
c = [];
p = [];
for idx=1:length(a)
c(:,idx) = generalized_normcdf(x,mu,a(idx),2);
p(:,idx) = generalized_normpdf(x,mu,a(idx),2);
end

clf
subplot(2,2,1)
plot(x,c);
title('cdf')

subplot(2,2,3)
plot(x,p);
hold on
set(gca,'colororderindex',1)
plot(mid(x),diff(c)./diff(x),'--');

for idx=1:length(b)
c(:,idx) = generalized_normcdf(x,mu,1,b(idx));
p(:,idx) = generalized_normpdf(x,mu,1,b(idx));
end
%c = generalized_normcdf(x,mu,1,b);
%p = generalized_normpdf(x,mu,1,b);
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
c = generalized_normcdf(x,par(1),par(2),par(3))+1e-3*randn(n,1);


% par = fit_generalized_normcdf(x,c);

x=2.2
c = generalized_normcdf(x,1,2,3);
generalized_norminv(c,1,2,3)


k = (1:1e4)';
n = 1e4;
%b = sqrt(2*log(k));
b1 = sqrt(2*log(k));
b = sqrt(2*log(k))-0.5*(log(4*pi*log(n)))./sqrt(2*log(n));

% cf. cramer
b2 = (sqrt(2*log(k)) - (log(log(n))+log(4*pi))/(2*sqrt(2*log(n)))); %/sqrt(2*log(n))
b3 = (sqrt(2*log(k)) - (log(log(n))+log(4*pi*log(2).^2))/(2*sqrt(2*log(n)))); %/sqrt(2*log(n))
norm(b-b2)

x = randn(k(end),n);
m = mean(cummax(x),2);
W = double(eulergamma);
%m(:,2) = b - W;
m(:,2) = b + W./b;
m(:,3) = b3 + W./b3;
m(:,4) = b;

subplot(2,2,1)
semilogx(m)
%[m,b1])
subplot(2,2,2)
rel=(m(:,2:end)-m(:,1))./m(:,1)
loglog(abs(rel))

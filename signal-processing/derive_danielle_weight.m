% 2016-09-18 10:43:33.450714597 +0200

syms n a
f = (n-1+2*a)^2 - n*(n-1+2*a^2)
a=solve(f,a)
a=simplify(a,'IgnoreAnalyticConstraints',true)
a=a(1)

af= @(n) -(2^(-1/2)*(n.^2 - n).^(1/2) - n + 1)./(n - 2)
a=af(1:1e4);
subplot(2,2,1)
plot(a)
subplot(2,2,2)
loglog(abs(a-(1-1/sqrt(2))))
a=af(10)
w=[a ones(1,9), a];
sum(w).^2/sum(w.^2)


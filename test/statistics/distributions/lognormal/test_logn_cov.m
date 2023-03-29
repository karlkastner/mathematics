% Wed 29 Mar 11:19:15 CEST 2023

n = 1e7;
p = 0.6;
mua = 0.2;
mub = 0.3;
sa = 0.5;
sb = 0.7;

za = randn(n,1);
zc = randn(n,1);
% correlate the variables
zb = (p*za + sqrt(1-p^2)*zc);

a = mua + sa*za;
b = mub + sb*zb;

corr(a,b)

ea = exp(a);
eb = exp(b);

% test covariance
'cov numeric'
cov(ea,eb)
mean(exp(a+b)) - mean(ea)*mean(eb)
'cov formula'
logn_cov(p,mua,mub,sa,sb)

% test correlation

corr(ea,eb)
logn_corr(p,mua,mub,sa,sb)

r = linspace(0,1)';
r(:,2) = logn_corr(r);
plot(r(:,1),r)

%if (0)
%mean(exp((1+p)*x).*exp(sqrt(1-p^2)*y_)), exp((1+p)^2/2+(1-p^2)/2)
%end

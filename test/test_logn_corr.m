% Fri 21 Jul 21:16:22 CEST 2023

mu = 0.2;
sd = 0.3;
mu2 = 0.5;
sd2 = 0.7;
r = 0.4;

[lmu,lsd]   = logn_moment2par(mu,sd);
[lmu2,lsd2] = logn_moment2par(mu2,sd2);


n = 1e7;
lx = randn(n,1);
ly = r*lx + sqrt(1-r^2)*randn(n,1);

corr(lx,ly)

%x = lognrnd(lmu,lsd,1e3,1);
%y = lognrnd(lmu2,lsd2,1e3,1);
x = exp(lmu + lsd*lx);
y = exp(lmu2 + lsd2*ly);
corr(x,y)
logn_corr(r,lmu,lmu2,lsd,lsd2)

% test mean
mean(x)
mean(y)
mean((x+y)/2)


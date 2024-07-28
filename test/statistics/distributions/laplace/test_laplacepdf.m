
x = innerspace(0,1e3,1e6)';
fx = fourier_axis(x);
df=fx(2)-fx(1);

mu = 1;
s  = 0.1;


S      = laplacepdf(fx,mu,s);
S(:,2) = laplacepdf(fx,mu,s) + laplacepdf(fx,-mu,s);

R = laplaceacf(x,mu,s);
R(:,2:3) = fft(S);
R(:,2:3) = R(:,2:3)./R(1,2:3);

w = laplacepdf_width(mu,s);

subplot(2,2,1)
plot(fx,S);
hline(0.5*max(S(:,1)))
vline(mu+w/2);

subplot(2,2,2)
plot(x,R)

% /home/pia/phd/src/lib/mathematics/statistics/distributions/laplace/laplacepdf_mean.m

h    = laplacepdf_entropy(mu,s);
fdx = (fx>=0);
lnS = log(S);
lnS(S==0) = 0;
h(2:3) = -sum(S(fdx,:).*lnS(fdx,:)*df);
h(1:3)

sd = laplacepdf_std(mu,s);
s2 = sum((fx(fdx)-mu).^2.*S(fdx,:).*df);
sd(2:3) = sqrt(s2);
sd(1:2)

ku = laplacepdf_kurtosis(mu,s);
s4 = sum((fx(fdx)-mu).^4.*S(fdx,:).*df);
ku(2:3) = s4./sd(1).^4;
ku

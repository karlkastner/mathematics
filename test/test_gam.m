% 2021-09-28 21:22:29.774643029 +0200 ../test_gam.m
% absolute
mu=2; sd=3; x=gamrnd(mu^2/sd^2,sd^2/mu,1e6,1); [mean(x),std(x)]

% relative
mu=2; sd=2; x=gamrnd(1/sd^2,sd^2*mu,1e6,1); [mean(x),std(x)/mu]
mu=2; sd=3; x=mu*gamrnd(1/sd^2,sd^2,1e6,1); [mean(x),std(x)/mu]

m=2; s=3; [a,b] = gam_moment2param(m,s); x = gamrnd(a,b,1e6,1); [mean(x),std(x)]
m = 1; s = 0.25; x = linspace(0,40,1e4)';  [a,b] = logn_moment2param(m,s); p = lognpdf(x,a,b); y = lognrnd(a,b,1e6,1); [a,b] = gam_moment2param(m,s); y(:,2) = gamrnd(a,b,1e6,1); p(:,2) = gampdf(x,a,b); y(:,3) = normrnd(m,s,1e6,1); p(:,3) = normpdf(x,m,s); plot(x,p),
[mean(y); std(y); skewness(y); kurtosis(y)] 

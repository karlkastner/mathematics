n = 1e6;
p = 0.1;
mu1 = 2;
mu2 = 3;
sd1 = 1./5;
sd2 = 1./7;

x = [normrnd(mu1,sd1,p*n,1); normrnd(mu2,sd2,(1-p)*n,1)];
mu = mean(x)
sd = std(x)
sk = skewness(x)

binormpdf_mean(p,mu1,mu2,sd1,sd2)
binormpdf_std(p,mu1,mu2,sd1,sd2)
sk_ = binormpdf_skewness(p,mu1,mu2,sd1,sd2)


x = linspace(-4,4)';

mu = 0;
sd = 1;

sk = [0,-1,-2,-4];
ku = [0,10,20,30];

clf
for idx=1:4
for jdx=1:4
y(:,1) = skew_generalized_normpdf(x,mu,sd,sk(idx),ku(jdx))
%y(:,2) = skew_generalized_normpdf(x,mu,sd,sk,25)
%y(:,3) = skew_generalized_normpdf(x,mu,sd,sk,50)
%y(:,4) = skew_generalized_normpdf(x,mu,sd,sk,100)
end
subplot(2,2,idx)
plot(x,y)

H = Histogram(y',inner2outer(x));
[H.mean,H.std,H.skewness,H.kurtosis]

end

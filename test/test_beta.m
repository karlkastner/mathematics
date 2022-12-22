% 2021-09-29 20:20:08.892427640 +0200
a = 2;
 k=[1,10,100,1e3,1e4,1e5], [chi2_mean(a), chi2_std(a), chi2_skew(a), chi2_kurt(a)], [k.*beta_mean(1,k-1);
 k.*beta_std(1,k-1);
 beta_skew(1,k-1);
 beta_kurt(1,k-1)], k=1e1;
 x=[1/2*chi2rnd(2,1e6,1),k*betarnd(1,k-1,1e6,1)];
 [mean(x);
std(x);
skewness(x);
 kurtosis(x)]

% 2016-02-29 13:36:34.401870431 +0100
% Karl Kastner, Berlin
n=100; p=(1:n)/(n+1);
sd = se_sample_quantiles(p);
plot(p,sd)



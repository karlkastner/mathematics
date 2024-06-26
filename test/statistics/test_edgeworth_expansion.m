% Thu 17 Aug 12:45:31 CEST 2023
% Karl Kastner, Berlin
% note that the lognormal distribution is not defined by its moments!


lc = 1;
Sc = 0.5;

[a,b] = logn_mode2par(lc,Sc);
%[a,b] = logn_moment2par(1,1);

mu = logn_mean(a,b);
sd = logn_std(a,b);
sk = logn_skewness(a,b);
k  = logn_kurt(a,b);

S = lognpdf(x,a,b);
S(:,2) = edgeworth_expansion(x,mu,sd);
S(:,3) = edgeworth_expansion(x,mu,sd,sk);
S(:,4) = edgeworth_expansion(x,mu,sd,sk,k);

plot(x,S)

% gaussian mixture:
% p, mu, sd, mu2, sd2


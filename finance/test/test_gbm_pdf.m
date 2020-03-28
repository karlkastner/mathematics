% Mon 20 Jan 01:27:20 +08 2020
r = 0.2;
s = 0.03;
T  = 5;
S0 = 10;
K  = 10.1;


m1 = quad(@(S_T) S_T.*gbm_pdf(T,S_T,r,s,S0),sqrt(eps)*K,K*10)
m2 = quad(@(S_T) S_T.^2.*gbm_pdf(T,S_T,r,s,S0),sqrt(eps)*K,K*10)
one = quad(@(S_T) gbm_pdf(T,S_T,r,s,S0),sqrt(eps)*K,K*10)
sd = sqrt(m2 - m1^2)

m1(2) = gbm_mean(T,r,s,S0)
sd(2) = gbm_std(T,r,s,S0)


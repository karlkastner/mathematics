% Fri 24 Feb 11:15:40 CET 2023

t = (0:0.01:100)/10;
S0 = 2;
r  = 0;
sigma = 1;

n = 10000;

S_sd = gbm_std(t(end),r,sigma,S0);

S    = gbm_simulate(t,r,sigma,S0,n);

S_sd
sqrt(mean(var(S(end,:))))

'geo'
S_sd = gbm_geostd(t(end),r,sigma,S0)
geostd(S(end,:))

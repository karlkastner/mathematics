% Fri 24 Feb 11:15:40 CET 2023

t = (0:0.01:100)/10;
S0 = 1;
r  = 1;
sigma = 0;

n = 10000;

S_mu = gbm_mean(t(end),r,sigma,S0)

S    = gbm_simulate(t,r,sigma,S0,n);

mean(S(end,:))

'geo'
S_mu = gbm_geomean(t(end),r,sigma,S0)
geomean(S(end,:))


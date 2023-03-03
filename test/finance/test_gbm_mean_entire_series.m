% Fri 24 Feb 11:15:40 CET 2023

t = (0:0.1:100)/10;
S0 = 1.1;
r  = 0.2;
sigma = 0.7;

%n = 1000000;

S_mu  = gbm_mean_entire_series(S0,r,sigma,t(end))

S_mu_ = gbm_mean(t,r,sigma,S0);
S_mu_ = mean(S_mu_)

%S    = gbm_simulate(t,r,sigma,S0,n);
%mean(S(end,:))


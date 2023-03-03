% Fri 24 Feb 11:15:40 CET 2023

t = (0:0.1:100)/10;
S0 = 1.1;
r  = 0.2;
sigma = 0.7;

%n = 1000000;

sd_S  = gbm_std_entire_series(S0,r,sigma,t(end))

sd_S_ = gbm_std(t,r,sigma,S0);
sd_S_ = sqrt(mean(sd_S_.^2))

%S    = gbm_simulate(t,r,sigma,S0,n);
%mean(S(end,:))


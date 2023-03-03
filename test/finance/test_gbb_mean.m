% Fri 24 Feb 11:15:40 CET 2023

t = (0:1:100)'/10;
S0 = 1.1;
r  = 0.5;
sigma = 1.7;

n = 100000;

S_mu = gbb_mean(S0,T,t,r,sigma);

S    = gbb_simulate(t,r,sigma,S0,n);
S_mu_ = geomean(S,2);

plot(t,[S_mu,S_mu_])

%mean(S(end,:))


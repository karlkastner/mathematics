

S0 = 2;
mu_S = 3
sd_S = 5
T = 7;

[r,sigma] = gbm_moment2par_entire_series(S0,mu_S,sd_S,T)

mu_S = gbm_mean_entire_series(S0,r,sigma,T)
sd_S = gbm_std_entire_series(S0,r,sigma,T)


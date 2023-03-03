mu_s = 2
sd_s = 3
S0   = 1.1
t    = 7



[r,sigma] = gbm_moment2par(mu_s,sd_s,S0,t)
mu_s = gbm_mean(t,r,sigma,S0)
sd_s = gbm_std(t,r,sigma,S0)



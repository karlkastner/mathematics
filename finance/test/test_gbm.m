% 2020-01-12 14:46:37.131732885

S0 = 1e3;
t  = 365;
n  = 10000;

% for one year
r     = 0.10
sigma = 0.2

% for one day
[r, sigma] = gbm_transform_time_step(r,sigma,1/365,1)
%r     = (1+r)^(1/365)-1
%sigma = sigma/sqrt(365)


S_me = gbm_median(t,r,sigma,S0) 
S_me_ = gbm_inv(t,0.5,r,sigma,S0)
S_mu = gbm_mean(t,r,sigma,S0)
S_sd = gbm_std(t,r,sigma,S0)

s = 100000;

f   = gbm_pdf(t,S_mu,r,sigma,S0)
q_  = gbm_cdf(t,S_me,r,sigma,S0)

pause

dt = 0.1;
t = (1:dt:365)';
%t = (1:(10*365))';
%t = (1:(1*365))';
S = gbm_simulate(t,r,sigma,S0,n);

S_mu(2) = mean(S(end,:))
S_me(2) = median(S(end,:))
S_sd(2) = std(S(end,:))

mean(min(S))
median(min(S))

r_ = [];
sigma_ = [];
S0_ =[];
for idx=1:n
	[r_(idx,1),sigma_(idx,1),S0_(idx,1)] = gbm_fit(t,S(:,idx));
%	[r_(idx,1),sigma_(idx,1),S0_(idx,1)] = gbm_fit(t,S(:,idx),'start');
%	[r_(idx,2),sigma_(idx,2),S0_(idx,2)] = gbm_fit(t,S(:,idx),'none');
%	[r_(idx,3),sigma_(idx,3),S0_(idx,3)] = gbm_fit(t,S(:,idx),'end');
%	[r_(idx,4),sigma_(idx,4),S0_(idx,4)] = gbm_fit(t,S(:,idx),'other');
end
[r__, sigma__, S0__] = gbm_fit(t,mean(S,2));
r = [r,mean(r_),r__]
serr(r_)./mean(r_)
sigma = [sigma,mean(sigma_),sigma__]
mu = r + 1/2*sigma.^2
S0 = [S0,mean(S0_),S0__]
%[r(2),sigma(2),S0(2)] = gbm_fit(t,S(:,1),'end')
%[r(2),sigma(2),S0(2)] = gbm_fit(t,S(:,1))


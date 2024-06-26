% Wed 15 May 11:02:45 CEST 2024
% Karl Kastner, Berlin
% maybe subtract mean and devide by std to improve scaling?
%
% - predict transition with inverse cdf
% - influence of rate of change
% - influence of noise level and spectrum 

m = 6;

rng(1)
p1 = rand();
p2 = 1-p1;
mu1 = -2+randn()
mu2 =  2+randn()
sd1 = lognrnd(1,1)
sd2 = lognrnd(1,1)
par = [p1;mu1;mu2;sd1;sd2];

n = 1e4;
x1 = normrnd(mu1,sd1,n,1);
x2 = normrnd(mu2,sd2,n,1);
n1 = round(n*p1);

% not: x = p1*x1 + p2*x2 !!!
x = x1; %x(1:n1) = x1(1:n1);
x(n1+1:end) = x2(n1+1:end);

% relation 3
[p1*moments(x1,m,0)+p2*moments(x2,m,0),moments(x,m,0)]


%mu = mean(par);
%[mu_(1),m
%par(1:2)
%mean(x)
%std(x)


[par_5,resn,res,exitflag,out] = g2_fit(x,[],5)
[par_6,resn,res,exitflag,out] = g2_fit(x,[],6)
%par_ = fit_bimodal(x,par_)
[par, par_5,par_6]
hist(x)

mu = [g2_moments(par,m),g2_moments(par_,m),moments(x,m,true)]


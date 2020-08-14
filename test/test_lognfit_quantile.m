p = [0.5,0.9];

x = [200,500];

param = lognfit_quantile(p,x)

x_ = logninv(p,param(1),param(2))
x_ = logninv(p,mean(log(x)),std(log(x)))



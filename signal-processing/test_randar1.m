% Thu 17 Jan 11:55:48 CET 2019
sigma = 1;
p = 0.9;
n = 20;
m = 1e6;
mu = 0;
x = randar1(sigma,p,n,m,mu);

x = std(x,[],2);

plot(x)


% 2021-06-30 22:07:59.509934832 +0200

% achieve same sum 1 over "n" time steps and sd of s

n = 2;
m = 1e4;
s = 0.1;
r = gamrnd(1/(s^2*n),s^2,1e4,n);
g=sum(r,2);
[mean(g),std(g)]

n=100;
g=gamrnd(1/(s^2*n),s^2,n,1e4);
mean(sum(g,2))
std(sum(g,2)), mean(std(g))

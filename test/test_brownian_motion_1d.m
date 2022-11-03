
n = 1e1;
m = 1e6;

y = bm_1d_fourier([n,m]);
y = y-y(1,:);
s = std(y,[],2);
x = (0:n(1)-1)'/n(1);
plot(x,[s,sqrt(x)])

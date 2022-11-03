% 2015-08-03 09:45:30.617567876 +0200
% Karl Kastner, Berlin


%xcorr

n     = 10^5;
rho1  = 0.25;
rho2  = 0.5;
rho12 = 0.75;
s1 = 2;
s2 = 3;
[x y]=randar1_dual(s1,s2,rho1,rho2,rho12,n);

%std(x)
%std(y)
%acf_man(x,2)
%acf_man(y,2)
%corr(x,y)
nc = 10;
C = xcorr(x,y,10);
size(C)

clf
N = (-nc:nc)';
plot(N,C/n,'.-')
mc = max(C)/n;
hold on

plot(N, mc*[rho2.^max(-N,0) rho1.^max(N,0)],'.r')

%c_ = rho12*sqrt(1-rho1^2)*sqrt(1-rho2^2)/(1-rho1*rho2)
%corr(x,y)
corr([x,y])
cov(x,y)


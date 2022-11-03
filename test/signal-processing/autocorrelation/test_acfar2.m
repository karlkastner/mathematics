% 2019-06-27 18:24:45.153535546 +0200
n = 1e6;
m = 1000;
r1=0.45; r2=0.45; x  = randar2(r1,r2,n); ro=roots([1 -r1 -r2])
A=[autocorr(x,m) acfar2(r1,r2,m+1) ro(1).^(1:m+1)']; % acfar1(ro(1),1e6,m) ];
subplot(2,1,1)
plot(A)
subplot(2,1,2)

plot( (A(:,2)-A(:,3))./(A(:,2)))



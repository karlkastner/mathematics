% 18:05:43.402621801 +0200
rho1 = 0.3;
rho2 = 0.5;
n = 10;

I = repmat((1:n)',1,n);
D = I-I';

s1 = sum(sum(rho1.^(-D).*(D < 0) + rho2.^D.*(D>0)))
s2 = 0.5*sum(sum(rho1.^abs(D) + rho2.^abs(D)))-n


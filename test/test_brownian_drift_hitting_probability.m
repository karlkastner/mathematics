% c.f. Abundo
% -> this still does not work when the threshold is smaller than the first value
m = 1e4;
dx = 1;
L = 1000;

s = 1;
s = 1;


b =1*1*2*pi/100;
x = (0:dx:L)';

% lower bound
a = 2;
% start value
y0 = 1;

% simulate
e = randn(length(x),m);
e(1,:) = 0;
y = y0 + b*x + s*sqrt(dx)*cumsum(e);

P1 = mean(y<=a,2);

P_ = mean(cummax(y)>a,2);
P = mean(cummin(y)>a,2);
P(:,2) = brownian_drift_hitting_probability(x,a,y0,b,s);

clf();
subplot(2,2,1);
% probability, that it ends below a
P1(:,2) = normcdf(a,b*x+y0,s*sqrt(x))
P1(:,3) = normcdf((a-y0-b*x)./(s*sqrt(x)))
plot(x,P1)
title('P(y(x))<a')

P2      = mean((cummin(y)<=a) & (y>=a),2);
P2(:,2) = exp(-y0*b);
P2(:,3) = exp(-2*y0*b)*normcdf(-(a-y0-b*x)./(s*sqrt(x)));
%normcdf(((b*x+y0) - a)./(s*sqrt(x))) 
subplot(2,2,2)
plot(P2)
title('P(y(X)>a), but P(y(x)<a,x<X)')


subplot(2,2,3)
plot(x,[mean(y,2),mean(cummin(y),2)])

subplot(2,2,4)
plot(x,P(:,1))
hold on
plot(x,P(:,2))
plot(x,P_)

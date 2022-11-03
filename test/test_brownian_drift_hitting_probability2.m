% Mon  4 Jul 15:55:25 CEST 2022
% Note: there seems to be an error in albundo et al
% -> sometimes probability large 1
% -> only valid for certain cases
a = 0;
b = 1;
t = 1;
s = 1;

%a = linspace(-10,10)';
a = 2;
y0 = 0;

t = linspace(0,10)';

m = 1e3;
% simulate
x = t;
e = randn(length(x),m);
e(1,:) = 0;
y = y0 + b*x + s*sqrt(dx)*cumsum(e);
P = p_passage(t,a,-b,y0);
P(:,2) = 1-mean(cummax(y)<a,2);

clf();
plot(t,P);
%plot(a,P)

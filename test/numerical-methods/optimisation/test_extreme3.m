% 2023-03-17 18:39:36.817347296 +0100
% x^2, x, 1
c = [ 3,2,1]
x = [-1,0,1]
y = polyval(c,x)

[y0,x0] = extreme3(x,y,2)

x = linspace(-1,1,1e3);
y = polyval(c,x);
clf
plot(x,y);
hold on
plot(x0,y0,'*')



% 2019-06-26 19:17:23.632595543 +0200
% Test
n=-1/4;
x = 2 + 3i; x.^(1/n)
m=max(abs(n),abs(1/n));
abs(x).^(1/n)*exp(1i*angle(x)/n+2i*pi*(0:abs(m)-1)'/m)

imroots(x)

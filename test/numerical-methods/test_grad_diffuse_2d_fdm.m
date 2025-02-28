% Wed  6 Nov 11:15:37 CET 2024
n  = 100*[1,1];
L  = 9*[1,1];
dt = 1;
e  = 1;
x  = linspace(-L(1)/2,L(1)/2,n(1))';
r = hypot(x,x');
y  = flat(sinc(r));
%y = flat(y+y');
delta_e = 0.1;
% y_ = A \ y;

[dy_de,D2] = grad_diffuse_2d_fdm(y,dt,e,n,L);

I = speye(prod(n));
dy_de(:,2)      = ( (I-dt*(e+delta_e)*D2)\y - (I-dt*(e-delta_e)*D2) \y)/(2*delta_e);

subplot(2,2,1)
imagesc(reshape(dy_de(:,1),n));
colorbar
subplot(2,2,2)
imagesc(reshape(dy_de(:,2),n));
colorbar
subplot(2,2,3)
imagesc(reshape(dy_de(:,1)-dy_de(:,2),n));
colorbar
%plot(x,[dy_de]);


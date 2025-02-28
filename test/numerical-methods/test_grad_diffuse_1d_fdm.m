% Wed  6 Nov 11:15:37 CET 2024
n = 100;
L = 9;
dt = 1;
e = 1;
x = linspace(-L/2,L/2,n)';
y = sinc(x);
delta_e = 0.1;
% y_ = A \ y;

[dy_de,D2] = grad_diffuse_1d_fdm(y,dt,e,n,L);

I = speye(n);
dy_de(:,2)      = ( (I-dt*(e+delta_e)*D2)\y - (I-dt*(e-delta_e)*D2) \y)/(2*delta_e);

plot(x,[dy_de]);


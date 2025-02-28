% Wed  6 Nov 11:15:37 CET 2024
n = 100;
L = 15;
dx = L/n;
%L = 9;
dt = 1;
e = 1;
x = linspace(-L/2,L/2,n)';
y = sinc(x);
delta_e = 0.1;
% y_ = A \ y;

[dy_de] = grad_diffuse_1d_spectral(y,dt,e,n,L);
[dy_de(:,2),D2] = grad_diffuse_1d_fdm(y,dt,e,n,L);
delta_e = 0.01;
[dy_de(:,3)] = (step_advect_diffuse_spectral(y,dt,0,e+delta_e,n,L) ...
		- step_advect_diffuse_spectral(y,dt,0,e-delta_e,n,L)) ./ (2*delta_e);

a = 0;
[dy_de(:,4)] =	( step_advect_diffuse_euler_implicit(dt,dx,n,y,a,e+delta_e) ...
		  - step_advect_diffuse_euler_implicit(dt,dx,n,y,a,e-delta_e) ) / (2*delta_e);
%I = speye(n);
%dy_de(:,2)      = ( (I-dt*(e+delta_e)*D2)\y - (I-dt*(e-delta_e)*D2) \y)/(2*delta_e);
for idx=1:size(dy_de,2)
subplot(2,2,idx)
plot(x,dy_de(:,idx));
end
%subplot(2,2,1)
%plot(x,dy_de(:,1));
%subplot(2,2,1)
%subplot(2,2,1)


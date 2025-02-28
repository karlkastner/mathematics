% Wed  6 Nov 11:15:37 CET 2024
n  = 100*[1,1];
L  = 9*[1,1];
dt = 1;
e  = 1;
x  = linspace(-L(1)/2,L(1)/2,n(1))';
r = hypot(x,x');
y  = flat(sinc(r));
%y = flat(y+y');
delta_e = sqrt(eps)
% y_ = A \ y;
dx = L./n;


ynext   = step_advect_diffuse_spectral(y,dt,[0,0],[e,e],n,L);
[dy_de] = grad_diffuse_2d_spectral(y,dt,e,n,L);

dy_de(:,2) = (   step_advect_diffuse_spectral(y,dt,[0,0],[e,e]+delta_e,n,L) ...
	       - step_advect_diffuse_spectral(y,dt,[0,0],[e,e]-delta_e,n,L) ) / (2*delta_e);

dy_de(:,3) = grad_diffuse_2d_fdm(y,dt,e,n,L);

a = [0,0];
[dy_de(:,4)] =	( step_advect_diffuse_euler_implicit(dt,dx,n,y,a,e*[1,1]+delta_e) ...
		  - step_advect_diffuse_euler_implicit(dt,dx,n,y,a,e*[1,1]-delta_e) ) / (2*delta_e);
%I = speye(prod(n));
%dy_de(:,2)      = ( (I-dt*(e+delta_e)*D2)\y - (I-dt*(e-delta_e)*D2) \y)/(2*delta_e);
dy_de = real(dy_de);

for idx=1:2
subplot(2,3,1+3*(idx-1))
imagesc(reshape(dy_de(:,1+2*(idx-1)),n));
axis equal
colorbar

subplot(2,3,2+3*(idx-1))
imagesc(reshape(dy_de(:,2+2*(idx-1)),n));
axis equal
colorbar

subplot(2,3,3+3*(idx-1))
imagesc(reshape(dy_de(:,1+2*(idx-1))-dy_de(:,2+2*(idx-1)),n));
axis equal
colorbar
end

%colorbar
%subplot(2,3,2)
%imagesc(reshape(dy_de(:,2),n));
%colorbar
%subplot(2,3,3)
%imagesc(reshape(dy_de(:,3),n));
%%imagesc(reshape(dy_de(:,1)-dy_de(:,2),n));
%colorbar
%%plot(x,[dy_de]);
%
%subplot(2,3,4)
%
if (0)
subplot(2,3,4)
imagesc(reshape(y,n))
subplot(2,3,5)
imagesc(reshape(ynext,n))
axis equal
subplot(2,3,6)
imagesc(reshape(dy_de(:,1)./dy_de(:,2),n));
colorbar
end

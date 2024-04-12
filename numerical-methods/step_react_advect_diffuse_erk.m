% Wed 13 Mar 13:34:35 CET 2024
% adopted from Kassam & Trefethen
% generalized to allow arbirary diffusion - advection
% an systems of multiple variables
% z : nx * ny * nvar
function [z,out] = step_react_diffuse_erk(dt,dx,a,e,nfun,z,isreal_)

% TODO no magic numbers
M  = 16;

n = size(z);
nx = n(1:2);
nvar = size(z,3);

%nx = obj.nx;
%nvar = obj.nvar;
%nvar = n(3);

% TODO, this can be precomputed
L_ = nx.*dx;
[fx, fy] = fourier_axis_2d(L_,nx);
fy = fy';

fmax =  0.5*nx./L_;
p = 2/3;
ax = abs(fx)>=fmax(1)*(p);
ay = abs(fy)>=fmax(2)*(p);

%fz = fft2(z0);
L = zeros(nx(1),nx(2),nvar);
for idx=1:nvar
% half linear step
% note: sign seems correct
L(:,:,idx) = (  ((-4*pi*pi*e(1,idx))*fx + 2i*pi*a(1,idx)).*fx ...
                  + ((-4*pi*pi*e(2,idx))*fy + 2i*pi*a(2,idx)).*fy ...
              );
end
E2 = exp(0.5*dt*L);
% full linear step
E = E2.*E2; 

% no. of points for complex means
r  = exp(1i*pi*((1:M)-0.5)/M);
% roots of unity
L  = L(:);
LR = dt*L(:,ones(M,1)) + r(ones(nx(1)*nx(2)*nvar,1),:);
Q  = dt*real(mean((exp(LR/2)-1)./LR,2));
f1 = dt*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
%f2 = dt*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3,2)); this was wrong, why?
f2  =dt*real(mean((4+2*LR+exp(LR).*(-4+2*LR))./LR.^3 ,2));
f3 = dt*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
%f1=h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
%f3=h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
f1 = reshape(f1,[nx(1),nx(2),nvar]);
f2 = reshape(f2,[nx(1),nx(2),nvar]);
f3 = reshape(f3,[nx(1),nx(2),nvar]);
L  = reshape( L,[nx(1),nx(2),nvar]);
Q  = reshape( Q,[nx(1),nx(2),nvar]);
% clear LR

	% step in time
	% TODO this can be reused from the previous step
	v = fft2(reshape(z,[nx(1),nx(2),nvar]));
	%**Nonlinear terms are evaluated in physical space**
	% Nonlinear evaluation. g(u,*)
	Nv = (nfun(flat((ifft2(v)))));
	Nv = reshape(Nv,[nx(1),nx(2),nvar]);
	Nv = fft2(Nv);
	%Coefficient 'a' in ETDRK formula
	a  = E2.*v + Q.*Nv;
	%Nonlinear evaluation. g(a,*)
	Na = (nfun(flat((ifft2(a)))));
	Na = reshape(Na,[nx(1),nx(2),nvar]);
	Na = fft2(Na);
	%Coefficient 'b' in ETDRK formula
	b  = E2.*v + Q.*Na;
	%Nonlinear evaluation. g(b,*)
	Nb = (nfun(flat(ifft2(b))));
	Nb = reshape(Nb,[nx(1),nx(2),nvar]);
	Nb = fft2(Nb);
	%Coefficient 'c' in ETDRK formula
	c  = E2.*a + Q.*(2*Nb-Nv);
	% Nonlinear evaluation. g(c,*)
	Nc = (nfun(flat((ifft2(c)))));
	Nc = reshape(Nc,[nx(1),nx(2),nvar]);
	Nc = fft2(Nc);
	% update
	v  = E.*v + Nv.*f1 + (Na+Nb).*f2 + Nc.*f3;
%out = struct();
%	out.va = v;
	% anti alias
	v(ax,:) = 0;
	v(:,ay) = 0;

z = ifft2(v);
if (isreal_)
	z=real(z);
end
z = flat(z);
%out.f1 = f1;
%out.f2 = f2;
%out.f3 = f3;
%out.E2 = E2;
%out.Q= Q;
%out.a = a;
%out.b = b;
%out.c = c;
%out.Nv = Nv;
%out.Na = Na;
%out.Nb = Nb;
%out.Nc = Nc;
%out.alias = ax | ay;
%rms(z)
%pause(1)
end % step_diffuse_erk


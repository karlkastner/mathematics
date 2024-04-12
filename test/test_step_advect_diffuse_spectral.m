% 2024-02-13 15:31:44.240868098 +0100

dt=1/400;
a = [0,0];
e = 0.1*[1,1];
%e(1) = 0;
%e(2) = 0;
L = 128*[1,1];
%L = 4*[1,1];
n = 2*L;
dx=L./n;
t=10;
z0=zeros(n);
z0(end/2,end/2)=1;
% [z,fg] = step_diffuse_spectral(t,dx,n,z0,e);

subplot(2,4,1)
imagesc(z0);

subplot(2,4,5)
cla
plot(z0(:,end/2));


z=z0;
z_ = z0;
zs = z0;
%nt = 1;
dt = t/10;
%t/nt;
nt = t/dt;
q = 0.5;
for idx=1:nt
	idx
	zs = step_advect_diffuse_spectral(dt,dx,a,e,zs);
	%z = step_advect_diffuse_implicit_q_fft(dt,dx,a,e,z,0.5);
	z  = step_advect_diffuse_implicit_q_fft(dt,dx,a,e,z,q);
	z_ = step_advect_diffuse_implicit_q(dt,dx,n,z_(:),a,e,q);
	%z = step_advect_diffuse_implicit_q_fft(dt,dx,a,e,z,0);
end
z_ = reshape(z_,n);
subplot(2,4,2)
imagesc(zs);

subplot(2,4,5)
hold on
plot(zs(:,end/2));



subplot(2,4,3)
imagesc(real(z));
subplot(2,4,4)
imagesc(real(z_));

subplot(2,4,5)
hold on
plot(z(:,end/2));
plot(z_(:,end/2));



%plot([z,z0]);


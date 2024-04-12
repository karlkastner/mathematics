%dt=1/400;
%t = 1;
dt = 0.01;
t = 10;

nt = t/dt;
L=1e3*[1,1];
n=L;
dx=L./n;
z0=zeros(n);
z0(1,1)=1;
z=z0;

ee = [0.1,100];
for edx=1:length(ee)
e=ee(edx)*[1,1];
zmin = 0;
for idx=1:nt;
	z = step_diffuse_spectral(dt,dx,n,z,e);
	zmin(idx,1) = min(z,[],'all');
end
zt = step_diffuse_spectral(t,dx,n,z0,e);
subplot(2,2,1+2*(edx-1))
plot([z(:,1),zt(:,1)]);
subplot(2,2,2+2*(edx-1));
plot((1:nt)*dt,zmin)
end
%plot(z(1,:))


% Thu 29 Feb 15:14:34 CET 2024

if (0)
n = 128;
z0 = zeros(n,1);
a = 1;
e = 0;
dt = 1;
z0(end/2) = 1;
q  = 1;
a = -0.9;
D1 = derivative_matrix_1_1d(n,n-1,-sign(a),'circular');
D2 = derivative_matrix_2_1d(n,n-1,2,'circular');
I = speye(n);
A = (I-(q*dt)*(a*D1 + e*D2));
z = [z0,z0];
m = 10;
for jdx=1:m
z(:,1) = step_diffuse_implicit_q_fft(dt,a,e,z(:,1),q);
%z(:,1) = NaN;
z(:,2) = A \ ((I+((1-q)*dt*(a*D1 + e*D2)))*z(:,2));
end

%D2=full(D2);
%[V,e]=eig(D2);
%F=fft(eye(n));
%e=diag(e);
%x=(2*pi*((0:n-1))/(n))';
% v,V, v'*V, o = 2*pi*[(0:n/2),1:n/2-1]; e = 2*(cos(o/n)-1) 
clf
plot([real(z),z0]);
rms(z(:,1)-z(:,2))
else

% 2D

for idx=1:3
switch (idx)
case {1}
 % e = [1,0];
 a = [1,0]; 
 e = [0,0];
case {2}
 %e = [0,1];
 a = [0,1]; 
 e = [0,0];
case {3}
 %e = [1/3,3];
 a = [3,3];
e = [0,0];
% e = [1/3,3];
% e = [1,1];
end

n = 16;
z0 = zeros(n,n);
dt = 0.25;
z0(end/2,end/2) = 1;
q = 0.5;
q = 1;

z = step_diffuse_implicit_q_fft(dt,a,e,z0,q);

[D1x,D1y,D2x,Dxy,D2y] = derivative_matrix_2d(n*[1,1],(n-1)*[1,1],2,{'circular','circular'});
[D1x,D1y] = derivative_matrix_2d(n*[1,1],(n-1)*[1,1],-1,{'circular','circular'});
I = speye(n^2);
D = a(1)*D1x + a(2)*D1y + e(1)*D2x+e(2)*D2y;
A = (I-(q*dt)*D);
z(:,:,2) = reshape(A \ ((I + ((1-q)*dt)*D)*z0(:)),n,n);

%D2=full(D2);
%[V,e]=eig(D2);
%F=fft(eye(n));
%e=diag(e);
%x=(2*pi*((0:n-1))/(n))';
% v,V, v'*V, o = 2*pi*[(0:n/2),1:n/2-1]; e = 2*(cos(o/n)-1) 
subplot(3,4,1+(idx-1)*4)
plot(real([z(:,end/2,1), z(:,end/2,2)]));
subplot(3,4,2+(idx-1)*4)
imagesc(real(z(:,:,1)))
colorbar
subplot(3,4,3+(idx-1)*4)
imagesc(real(z(:,:,2)))
colorbar
subplot(3,4,4+(idx-1)*4)
imagesc(real(z(:,:,1)-z(:,:,2)))
colorbar
rms(z(:,:,1)-z(:,:,2),'all')
end
end


% Mon 13 Jun 13:01:16 CEST 2022

n = [10,1]
n = [20,20];
n = [5,5];
f = 5;
m=1e3;
L = 1*[1,1];
y = zeros(n(1),n(2),m);
x=(0:n(1)-1)'/n(1);
tic
[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,{'dirichlet','dirichlet'});
D2 = D2x+D2y;
% for idx=1:m
if (0)
y = bm_2d_fourier([n,m]);
else
[y,C,sC] = brownian_motion_2d_cov([n,m]);
end
% end;
toc/m
k = n/2;
k = [1,1];
%k =n(1);
 y=y-y(k(1),k(2),:);
s=(std(y,[],3));

subplot(3,3,1);
s_=[s(:,k(2)),sqrt(abs(x-x(k(1))))];
%s_=[s(k(1),:)',s(:,k(2)),sqrt(abs(x-x(k(1))))];
%s_=[s(k,1),s(1,k)',diag(s),sqrt(abs(x-x(k)))];
sc=1./(1/x(end)*s_(end,:))
plot(x,s_.*sc)
subplot(3,3,2);
imagesc(s);
axis square
colorbar
%y = y-y(1,:);
% plot(x,[std(y,[],2),sqrt(x)])
subplot(3,3,3)
imagesc(mean(y,3))
caxis([-1,1])
colorbar
axis square

s = 0.1;
b = sin(2*pi*f*(x + s*y));
subplot(3,3,4)
imagesc(b(:,:,1))

S = mean(abs(fft(fft(b,[],1),[],2)).^2,3);
subplot(3,3,5)
plot(S(1,:))
subplot(3,3,6)
%plot(S(:,2))


s = [0.05,0.1,0.2];
for idx=1:3
	b = sin(2*pi*f*(x + s(idx)*y(:,:,1)));
	subplot(3,3,6+idx)
	imagesc(b)
end
colormap gray

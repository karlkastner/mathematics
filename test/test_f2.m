% Tue 14 Sep 16:22:37 CEST 2021

n = 1e4;
N = [1e2,1e3,1e4,1e5,1e6];
figure(1)
clf
k = 10;
max_S = [];
IS = [];
sb=0.5;
se=eps;
L=100;
for idx=1:length(N)
n_ = N(idx);
%n_ = n;

x = L*(0:n-1)'/n;
y = [sqrt(2)*sb*sin(2*pi*k*x),se*randn(n,1)];
%sy = std(y);
%f = sqrt(2*L)./(sy*n)*fft(y);
%S = f.*conj(f);
y(:,3) = sum(y(:,1:2),2);
S = periodogram(y,L,n_);
%S(:,3) = mean(S(:,1:2),2);
df_=1/L;
max_S(idx,:) = max(S)
figure(1)
subplot(2,3,idx)
cla
if (n_ < n)
	%fx  = fourier_axis(x);
	x = linspace(0,n_/n*L,n_);
	fx = fourier_axis(x);
%	fdx = find(fx>=0);
%	fx  = fx(fx>=0);
	df(idx,1) = fx(2)-fx(1);
	%plot(fx(1:n_/2),S(end/2+1:end,:))
	plot(fx,S)
%(end/2+1:end,:))
else
	% padd zeros in space
	x = linspace(0,n_/n*L,n_);
	fx = fourier_axis(x);
	df(idx,1) = fx(2)-fx(1);
plot(fx,S)
end
IS(idx,:) = sum(S(n_/2+1:end,:))*df(idx);
vline(k,'color','b');
%fx = (0:n_/2-1)'*2*L/n_;
%df = fx(2)-fx(1);
l=length(S);
l_ = n; 
S(end/2+1:end,:) = 0;
R_ = L*real(ifft(S(:,[1,3])));
%if (l>n)
%end

%R_ = autocorr_fft(y,[],true);
%xlim([0,l/2])
R(idx,1) = R_(1);
figure(2)
subplot(2,3,idx)
cla
plot(R_);
title([n,n_])
xlim([0,l/2])
end
IS
m = 10;
Sb = bartlett_periodogram(y,L,m,n);
sum(Sb(end/2+1:end,:))/L


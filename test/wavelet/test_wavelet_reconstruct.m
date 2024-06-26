% 2016-07-09 15:59:14.978537949 +0200
% Karl Kastner, Berlin

m = 2000;
dt = 0.01;
t = dt*(1:m)';

f0 = [1];

n  = 7;
winstr = 'flattopwin';

% generate a sine wave
x = 0;
for idx=1:length(f0)
	x =     sin(2*pi*f0(idx)*t);
end

% forward transform
w = wavelet_transform(x,dt,f0,n,winstr);

% inverse transform
xr = wavelet_reconstruct(w,dt,f0,n,winstr);

clf
subplot(2,1,1)
plot(real(abs(w)))
subplot(2,1,2)
plot([x]); %, real(circshift(xr,[2 0]))],'.');
hold on
plot([xr],'.'); %, real(circshift(xr,[2 0]))],'.');


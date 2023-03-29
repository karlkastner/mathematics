% Thu 24 Jun 14:58:46 CEST 2021

i=1;
 nf=60;
 n = 4e2;
 x = randn(n);
y = bandpass2d_implicit(x,0.24,i);
% y = lowpass2d_implicit(x,0.2,i);
% y =x;
%y = x;
 f=fft2(y)/sqrt(numel(f));
 f=abs(f).^2;
 f=fftshift(f);
f=lowpass2d_implicit(f,0.1,i);
%f = meanfilt1(meanfilt1(f,10)',10)';
 imagesc(f);
 axis equal;
 var(x(:)), var(y(:)), rms(x(:)-y(:));
 colormap(parula(10))
colorbar

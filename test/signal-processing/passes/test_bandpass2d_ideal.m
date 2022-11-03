% Wed 27 Apr 17:19:16 CEST 2022
% This demonstrates that isotropic (spotted/gapped) patterns should not be
% analyzed along transects, bc the SD along the transect differs from the radial 2D SD


n=2*800+1;
 dx=0.5;
 L=2;
L_ = n*dx;
df = 1/L_;


x = (1:n)*dx;
fx = fourier_axis(x);

 p = 2;
 x = randn(n);
 [y,S,r]=bandpass2d_ideal(x,L,dx,p);
 imagesc(x);
 y_ = lowpass2d_ideal(x,L,dx,p);
 y_ = highpass2d_ideal(y_,L,dx,p);
 rms(flat(y-y_))
figure(1)
clf
subplot(2,2,1)
imagesc(y)
subplot(2,2,2)
imagesc(y_)

[S2b] = spectral_density_bandpass2d_ideal(size(x),dx,L,2);
S2h = spectral_density_highpass2d_ideal(size(x),dx,L,2);
S2l = spectral_density_lowpass2d_ideal(size(x),dx,L,2);
%S2 = fftshift(S2);

S2b = fftshift(S2b);

figure(2)
subplot(2,3,1)
imagesc(S2l)
subplot(2,3,2)
imagesc(S2h)
subplot(2,3,3)
imagesc(S2b)

figure(3)
subplot(2,2,3)
clf;
% nb: we cannot call periodogram1 here, as the columns have to be weighted by their spectral energy
%S1 = periodogram(y-0*mean(y),dx*n);
%Sx = fftshift(mean(S1,2));
Sx = mean(abs(fft(y)).^2,2);
Sx = Sx/(sum(Sx)*df);
%S1 = periodogram(y'-0*mean(y'),dx*n)';
%Sy = fftshift(mean(S1,1)');
%S_ = Sx+Sy;
S_ = 2*fftshift(Sx);

plot(fftshift(fx),S_); %/mean(S_))
hold on
S2b = S2b/(sum(S2b(:))*df^2);
S_ = S2b((n+1)/2,:)';
plot(fftshift(fx),[S_]); %/mean(S_));
%plot(fftshift(fx),[S_,sqrt(S_)]); %/mean(S_));

S = abs(fft2(y)).^2;
S = fftshift(S);
S = S/(sum(S(:))*df^2);
S = trifilt2(S,21);
S_ = S((n+1)/2,:)';
plot(fftshift(fx),[S_]); %/mean(S_));

R=ifft2(fftshift(S2b));
R = real(R);
 S=abs(fft(R)).^2;
 S = mean(S,2);
 S = sqrt(S);
 S=2*S/(sum(S)*df);
 hold on;
 plot(fx,real(S),'--')

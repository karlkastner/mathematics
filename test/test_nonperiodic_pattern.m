% an aperiodic (short) domain result in an underestimate of regularity (by about 20%)
% for bandpass
% for brownian phase: 
% padding the mirror image of the pattern increases the spectral resolution and improves the regularity estimate
% -> however, overestimates correlation for cross-patterns
 clf;
n=1e3;
 f0=0.1;
L = 10/f0;
df = 1/L;
 m=1e4;
p=10;
fx=fourier_axis(L,n);
S = spectral_density_bandpass_continuous(fx,f0,p);
%S = spectral_density_brownian_phase_across(fx,0.2);
% S = spectral_density_brownian_phase(fx,f0,0.3);
 b = real(ifft(sqrt(S).*randn(n,m)));
 S=mean(abs(fft(b)).^2,2);
 S=S/(sum(S)*df);
subplot(2,2,1)
 plot(fx/f0,S);
hold on
subplot(2,2,2)
plot(real(ifft(S)))
 b = b(1:end/2,:);
 S=mean(abs(fft(b)).^2,2);
 S=S/(2*sum(S)*df);
 hold on;
subplot(2,2,1)
 plot(fx(1:2:end)/f0,S);
subplot(2,2,2)
plot(real(ifft(S)))
w = trapwin(1:n)';
w = tukeywin(n);
w = 1;
fb = flipud(b);
 b = [b;
      fb];
%      fb(1:end/2,:)];
  S=mean(abs(fft(w.*b)).^2,2);
L_ = 1*L;
df_ = 1/L_;
n_ = 1*n;
 S=S/(sum(S)*df_);
 hold on;
subplot(2,2,1)
fx_ = fourier_axis(L_,n_);
 plot(fx_/f0,S);
%set(gca,'yscale','log');
xlim([0,10])
ylim([1e-7,25]);
%max(ylim)])
subplot(2,2,2)
plot(real(ifft(S)))

fx=fourier_axis(2*L,2*n);
%S = spectral_density_brownian_phase_across(fx,0.2);
%S = spectral_density_brownian_phase(fx,f0,0.3);
S = spectral_density_bandpass_continuous(fx,f0,p);
b = real(ifft(sqrt(S).*randn(2*n,m)));
b = b(1:end/2,:);
S=mean(abs(fft(b)).^2,2);
S=S/(sum(S)*df);
subplot(2,2,1)
plot(fx(1:2:end)/f0,S);
subplot(2,2,2)
plot(real(ifft(S)))

w = cvec(tukeywin(n,1));
S=mean(abs(fft(w.*b)).^2,2);
S=S/(sum(S)*df);
subplot(2,2,1)
plot(fx(1:2:end)/f0,S);
legend('L=8 lambda_c periodic','l = 4 lambda_c non-periodic','l = 4 lc mirrored','l = 8, non-periodic')
legend('L=8 lambda_c periodic','l = 4 lambda_c non-periodic','l = 4 lc mirrored','l = 8, non-periodic','l8 non-p windowed')
subplot(2,2,2)
plot(real(ifft(S)))

n = 1e5;
fx=fourier_axis(linspace(0,1,n));
c  = [1, -0.99, 0.6],
c  = [1, -1.1596,-0.3228]
a  = c(3); b = c(2);
fm = frequency_ar2(c)
r = roots(c)
(b+[1,-1]*sqrt(b^2 - 4*a))/2
%1/(2*pi)*acos(c(2)/(2*sqrt(c(3))))
y = randn(n,1)';
y=filter(1,c,y)';
%x_=x-mean(x);
subplot(2,2,1);
plot(y);
xlim([0,1000]);
subplot(2,2,2);
%S = periodogram_bartlett(y-mean(y),1,10,n);
S = periodogram(y-mean(y),1);
S(:,2) = spectral_density_ar2(fx,c);
S = S./sum(S);
[mv,mdx] = max(S);
plot(fx,S);
%vline(max(fx)*fm,'colo','r');
xlim([0,n/2]);
%plot([0.3*x/std(x),autocorr_fft(x)]) 
fx(mdx(2))
fx(mdx(2))/max(fx)
R = autocorr_fft(y);
subplot(2,2,3)
plot(R)
ar2_acf2c(R(2),R(3))

sys = ar(y,2);
[num, den] = tfdata(sys);
cc_ = den{1};
cc_ = cc_(2:3)'

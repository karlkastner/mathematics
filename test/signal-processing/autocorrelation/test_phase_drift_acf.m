% Wed  1 Dec 11:51:16 CET 2021

n = 1e5+1;
L = 1000;
fc = 1;
s  = 0.1;
nb = 8e2;
x  = linspace(-L/2,L/2,n)';
x  = circshift(x,-(n-1)/2);
fx = fourier_axis(x);

dx = L/n;
ex = s*sqrt(dx)*cumsum(randn(n,1));

y = cos(2*pi*fc*(x+ex)) + 0*sin(2*pi*fc*ex);

R = autocorr_fft(y);

R(:,2) = acf_brownian_phase(x,fc,s);
p = s;
dx = L/n;
p = pi*s^2/dx;
S = spectral_density_lorentzian(fx,fc,p);
R(:,3) = real(ifft(S));
R = R./max(R);
%(1,:);

figure(1);
clf

subplot(2,2,1)
plot(x,R);

fdx = round(1/fc/dx);
vline(x(fdx))

log(R(fdx,:))/x(fdx)


xlim([-10,10]/fc);
R(abs(x)>10)=0;
S = [fft(R(:,1:2)),S,periodogram_bartlett(y,L,nb,n)];
S = [S, spectral_density_brownian_motion(fx,fc,s)];
S = real(S);
S = S./sum(S);


subplot(2,2,2)
plot(fx,S);
legend('from R exp','from R analytic','Lorentzian','Bartlett','BM')

subplot(2,2,3)
plot(x,R./cos(2*pi*fc*x))
ylim([-0.5,1.5]*1);
xlim([-10,10])


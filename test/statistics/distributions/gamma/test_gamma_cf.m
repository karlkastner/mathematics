% 2023-12-04 13:29:05.710697445 +0100

L = 30;
n = L^2;
[a,b] = gamma_mode2par(1,1);
x = linspace(0,L,n)';
fx=fourier_axis(x);

pdf=gampdf(abs(fx).*(fx>0),a,b);
cf = gamma_cf(abs(fx),a,b);
R=fft(pdf);
R=R/R(1);

subplot(2,3,1);
plot(fx,pdf);

subplot(2,3,2);
plot(x,[real(R),real(cf)]);

subplot(2,3,3);
plot(x,[-imag(R),imag(cf)],'--')   


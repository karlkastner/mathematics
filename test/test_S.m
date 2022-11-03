% 2022-06-01 03:25:12.266894650 +0200

 L=10;
 df=1/L;
n = 1e3;
a = 1/2;
e = 0.5*0.35*brownian_noise_2d_fft([L/a,a*L],n);
std(e(:))
%e = trifilt2(e,7);
x = linspace(0,L,n)';
fx = fourier_axis(x);
b=(cos(2*pi*(e+x)));
b = (1+b).^2;
b = b-mean(b(:));

figure(1)
clf;
subplot(2,2,1)
e_ = e-e(:,1);
 plot(mean(e_.^2));
 e_ = e-e(1,:);
 hold on;
 plot(mean(e_.^2,2))
subplot(2,2,2)

a=autocorr_fft(b);
 plot(mean(a,2));
 hold on;
 a = autocorr_fft(b');
 plot(mean(a,2)) 


figure(2)
clf

% b = e;
 S = abs(fft2(b)).^2;
 S=(fftshift(S));

 subplot(3,2,1);
 imagesc(b) %>mean(b))
axis equal
axis square
 subplot(3,2,2);
 S = trifilt2(S,3); 
% imagesc(log10(S))
 imagesc(fftshift(fx),fftshift(fx),S);
axis(3*[-1,1,-1,1])

 subplot(3,2,3);
 x=linspace(-L/2,L/2,n)';
 x=x+-0.5*(x(2)-x(1));
 S_=cvec(mean(S));
 S_=S_/(sum(S_)*df);
 %fun = @(x,p) p(1)*exp(-p(2)*abs(x));
 %fun = @(x,p) p(1)*abs(x).^-p(2);
% fun = @(x,p) p(1)*p(2).^abs(x);
 fun = @(x,p) p(1)./(1+p(2)*abs(x));
%abs(x).^-p(2);
%exp(-p(2)*abs(x));
 p = lsqnonlin(@(p) fun(x,p)-S_, [1,1]);
% S_(:,2)=exp(-0.02*abs(x));
S_(:,2) = fun(x,p);
 S_=S_./sum(S_);
 plot(S_);
subplot(3,2,5)
plot(real(ifft(ifftshift(S_))));

 subplot(3,2,4);
 S_=(mean(S'));
  S_=S_/(sum(S_)*df);
 plot(S_)
subplot(3,2,6)
plot(real(ifft(ifftshift(S_))));

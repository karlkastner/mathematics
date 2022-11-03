% S = Sp + Ss, as independent

n  = 1e5;
m  = 1e1;

sb = 2;
se = 1;
L = 100;
f = 10;
L_w = 0.25;
L_w = 1;

x = linspace(-L/2,L/2,n)';
dx = x(2)-x(1);
fx = fourier_axis(x);
df = fx(2)-fx(1);
wfun = @(x) exp(-1/2*(x/L_w).^2);
w    = wfun(x);
%w  = tukeywin(n,1);
w  = w/mean(w);
yp = sb*w.*(1+cos(2*pi*f*x));
yp = yp/std(yp);

ye = se*w.*randn(n,m);
%ye = ye./std(ye);

%yp = yp-mean(yp);
ye = ye-mean(ye);

%mean(std(ye))
ye = se*ye./std(ye);
%mean(std(ye))

y = ye + yp;

nflag = true;
p_ = @(y) abs(fft(y)/n).^2;
if (nflag)
	p = @(y) periodogram(y,L,n);
else
	p = p_;
end

S   = mean(p(y),2);
Sp  = p(yp);
Se  = mean(p(ye),2);
% re-normalize

S_   = mean(p_(y),2);
Sp_  = p_(yp);
Se_  = mean(p_(ye),2);

s  = mean(sum(y.^2)*dx,2);
sp = sum(yp.^2)*dx;
se = mean(sum(ye.^2)*dx,2);
[s,sp,se]

%Sp_ = p_(yp);
%[sp sum(Sp_*df)]
%Se_ = mean(p_(ye),2);
%[se sum(Se_*df)]
%S_ = Sp_ + Se_;
%sum([Sp_,Se_,S_])*df
%rms(diff(S,[],2))./rms(S(:,1))

% only when normaliyed
%if (nflag)
se = se;
S(:,2) = 1/(se^2 + sp^2)*(sp^2*Sp + se^2*Se);
S_(:,2) = (Sp_ + Se_);

figure(1);
clf
subplot(2,2,1)
plot([S])
xlim([0,2*f]);

subplot(2,2,2)
%S__ = [mean(abs(fft(y)).^2,2),  (mean(abs(fft(yp)).^2 + abs(fft(ye)).^2 + 0*2*real(fft(ye).*conj(fft(yp))) + 0*conj(fft(ye)).*fft(yp),2))];
plot(S_)
xlim([0,2*f]);

subplot(2,2,3)
plot([Sp,Se])
%0.5*sp^2*Sp])
xlim([0,2*f]);
title('S_p')

subplot(2,2,4)
%plot([Se_,0.5*se^2*Se])
plot([Sp_,Se_])
xlim([0,2*f]);
title('S_p')

figure(2)
clf
subplot(2,2,1)
plot(yp)
R = autocorr_fft(yp);
R(:,2) = wfun((x-x(1))/sqrt(2)).*(2+1*cos(2*pi*f*(x-x(1))))/3;
subplot(2,2,2)
plot(R)

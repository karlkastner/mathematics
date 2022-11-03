% Wed  6 Oct 17:45:49 CEST 2021
% dx = [1,2,3,4]
 mm = 10;
 nn = 1e4;
 LL = [2.5,5,7.5,10,20];
% l = 1e2;

f1 = 10;
F2 = f1*[1,2,3,4,5];
%F2 = f1;
nf = 1;
fmin =0;
fmax = inf;

D =zeros(length(F2),length(LL));
p =zeros(length(F2),length(LL));
R = [];
l = 1;
for idx=1:l
for jdx = 1:length(LL)
 for idx=1:length(F2)
L = LL(jdx);
n = nn;
%n = L*nn/10;
%m = L*mm/10;
m = 10;
x = linspace(0,L,n)';
% xi = x;
% xi = linspace(0,L,10*n);
% f  = 10*(1+0*x);
 f2 = F2(idx);
 rho1 = filter_f0_to_rho(f1,L,n);
 rho2 = filter_f0_to_rho(f2,L,n);
 a = 1;
 %y = a.*sin(2*pi*f.*x) + 0*randn(n,1);
% y = a.*(x.*sin(2*pi*f.*x) + (1-x).*sin(2*pi*f2*x));
% y = y+randn(n,1);
 e  = randn(n,1);
 y1 = bandpass1d_implicit(e,rho1,1,true);
% y1 = highpass1d_implicit(e,rho1,1,true);
% y1 = lowpass1d_implicit(e,rho1,1,true);
% e_  = randn(n+2*nf,1);
% y1 = meanfilt1(e_,nf);
% y1 = y1(nf+1:n+nf);
 y2 = bandpass1d_implicit(e,rho2,1,true);
 y1 = y1/rms(y1);
 y2 = y2/rms(y2);
 %y(1:n/2) = y1(1:n/2);
 y = y1.*(x/L) + (1-x/L).*y2;
 fx=fourier_axis(x);
% y = interp1(x,y,xi,'linear');
 [p_,D_,pp,ratio,SS,rr] = periodogram_test_stationarity(y,m,L,fmin,fmax);
p(idx,jdx) = p(idx,jdx) + p_;
D(idx,jdx) = D(idx,jdx) + D_;
R(idx,jdx) = rms(pp(:,2)-pp(:,1));
%R(idx) = rms(ratio(:,2)-ratio(:,1));
end;
end
end
%plot(fx,periodogram_bartlett(y,1,m,n));
% xlim([0,50]);
% mean(p), mean(D), plot(SS);
%mean(p)
%mean(D)
%mean(R)
D = D/l;
p = p/l;

subplot(2,2,1)
plot(F2/f1,p)
legend(num2str(cvec(LL)))
hline([0.05])
subplot(2,2,2)
plot(F2/f1,D)
%plot(pp)

subplot(2,2,3)
cla
plot(fx,periodogram_bartlett(y(1:end/2),1/2,m,n));
hold on
plot(fx,periodogram_bartlett(y(end/2+1:end),1/2,m,n));

subplot(2,2,4)
plot(ratio)
%cla
%plot(y)
%plot(sort(D))
%plot(std(SS./mean(SS,2),[],2))
%S = abs(fft(y)).^2;
%S = S(1:end/2);
%SS = reshape(S,m,n/2/m)';
%plot(mean(SS,2))
%plot(std(SS./mean(SS,2),[],2))
%plot(mean(SS,2))
%plot([mean(rr,2),std(rr,[],2)])

% Wed  5 Mar 10:44:28 CET 2025

aa = [0.25,0.5,0.75];
L = 1e2;
n = 1e4;
%x=(0:n-1)'*L/n;
x = fourier_axis(n/L,n);
fx = fourier_axis(L,n);
fmax = max(fx);
clf
for idx=1:length(aa)
a=aa(idx);
S = abs(fx).^-a;
S = S.*(fx>0);
df = fx(2)-fx(1);
S(1)=0;	
S = S/sum(S*df);
R = real(fft(S));

subplot(2,3,idx)
loglog(fx,S)

%subplot(2,3,3+idx)
subplot(2,3,4)
%cla
loglog(x,R);
hold on
%w=real(expint(2i*x*pi*fmax));
w = 1;
plot(x, R(2)*w.*(x/x(2)).^-(1-a));


end

%R(:,2) = poweracf(x,a);

%subplot(2,2,1)
%plot(S)

%subplot(2,2,2)
%plot(R)

%R=R/R(1);
%R(:,2) = 1./(x.^(1-a));
%R=R./R(2,:);
%plot(R)

% 0.9	 0.06
% 0.5    0.217


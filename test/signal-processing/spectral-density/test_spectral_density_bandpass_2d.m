% Tue 23 Aug 13:57:17 CEST 2022

fx = linspace(0,1e5,1e6)';
p  = [1.01,4,16,64];
p = [1.5,3];
%p  = [1.1,2,4,8,16,32,64];
fc = [1,2];
S = [];
Sn = [];
k = 0;
Sc = [];
reg = [];
for idx=1:length(p)
for jdx=1:length(fc)
	k = k+1;
	S(:,k)  = spectral_density_bandpass_2d(fx,fc(jdx),p(idx),false);
	[Sn(:,k),Sc(k,1)] = spectral_density_bandpass_2d(fx,fc(jdx),p(idx),true);
	reg(k,1) = sqrt(Sc(k,1))*fc(jdx);
end
end

IS = cumsum(2*pi*mid(fx).*mid(Sn).*diff(fx));

figure(1);
clf
subplot(2,2,1)
%loglog(fx,S);
plot(fx,S);
xlim([0,10])
subplot(2,2,2)
plot(mid(fx),IS)
legend(num2str(reg));

subplot(2,2,3)
plot(fx,Sn)
xlim([0,10])

if (0)
% Wed 27 Apr 16:38:35 CEST 2022

Lf  = 2.5;
n = 4*[100,100]+1;
dx = 0.1;
L = dx*n(1);
x = linspace(-L/2,L/2,n(1));
fx = fourier_axis(x);
fr = hypot(fx,fx');
fc = 1/Lf;
df = 1/L;

p = 0.99;

[S,R,r] = spectral_density_bandpass2d_ideal(n,dx,Lf,p);
s = (df^2*sum(sum(S)))
S   = S/s;

S_ = spectral_density_bandpass_2d(fr,fc/(2*pi),p);

% test integral to 1
sum(sum(S_))*df^2

S  = fftshift(S);
S_ = fftshift(S_);

S   = S(:,(n(1)+1)/2);
S1_ = S_(:,(n(1)+1)/2);
%S   = S/max(abs(S));
%S1_ = S1_/max(abs(S1_));

plot(fftshift(fx),real([S,S1_]));
end

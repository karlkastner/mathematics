% Fri 22 Apr 13:27:51 CEST 2022

a = 1;
L = 10;
n = 999;
k = (n+1)/2;

x   = linspace(-L,L,n);
dx  = x(2)-x(1);
rho = exp(-a*dx)
% a = -log(rho)/dx

fx = fourier_axis(x);
r = hypot(x,x');
fr = hypot(fx,fx');
R_lp = exp(-a*r);
df  = (fx(2)-fx(1))^2;

% no square here
S_lp =((abs(fft2(R_lp))));

R_hp =  R_lp;
R_hp = -R_hp/sum(R_hp(:));
%f = fftshift(f);
R_hp(k,k)=R_hp(k,k)+1;
f   =R_hp/R_hp(k,k);
%f = ifftshift(f);
S_hp=((abs(fft2(f))));

S_bp = S_lp.*S_hp;
S_bp= S_bp*1/(sum(S_bp(:))*df^2);

R_bp = fftshift(real(ifft2(S_bp)));

S_lp = fftshift(S_lp);
S_hp = fftshift(S_hp);
S_bp = fftshift(S_bp);

clf
subplot(2,3,1);
plot(x,R_lp(k,:));

subplot(2,3,2);
plot(x,R_hp(k,:));

subplot(2,3,3);
plot(x,R_bp(k,:));


subplot(2,3,4);

S = cvec(S_lp(k,:)/sum(S_lp(k,:)));
normalize = @(x) x./sum(x);
fx = fftshift(fx);
Sfun = @(rho) normalize(spectral_density_lowpass(fx,rho,1,dx,'rho'));
%plot(fftshift(S/sum(S)))
%lrho = lsqnonlin( @(lrho) Sfun(exp(-rho*dx)) - S, -log(rho)/dx) 
%rho = exp(-lrho*dx)
rho0 = rho;
rho = lsqnonlin( @(rho) Sfun(rho) - S, rho)
%rho = exp(-lrho*dx)
%rho = 1.02*rho;
%rho = 1.03*rho;
%S(:,2) = Sfun(rho);
%S(:,3) = Sfun(1.01*rho);
%S(:,3) = Sfun(1.02*rho0);
%S(:,5) = Sfun(1.03*rho);
S(:,2) = spectral_density_lowpass_2d(fx,1/a,1);
S = S./sum(S);
%rms(S-S(:,1))

plot(fx,S)
%S_lp(k,:)/sum(S_lp(k,:)));

subplot(2,3,5);
S = cvec(S_hp(k,:));
S(:,2) = spectral_density_highpass_2d(fx,1/a,1);
S=S./sum(S);
plot(S);

subplot(2,3,6);
S = cvec(S_bp(k,:));
S(:,2) = spectral_density_bandpass_2d(fx,1/a,1);
S = S./(sum(S)*df^2);
plot(S);

p = 2;
fr = fftshift(fr);
S = spectral_density_bandpass_2d(fr,1/a,p);
Sc=sd_bp_2d_Sc(1/a,p)
sum(sum(S))*df^2/Sc
% test radial integral
S_ = S_bp(500:end,500); fr_ = (0:length(S_)-1)'*df; sum(2*pi*S_.*fr_)*df

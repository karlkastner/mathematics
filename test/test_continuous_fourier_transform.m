% Wed 15 Sep 13:26:31 CEST 2021

n  = 1e4+1;
nf = n; 

se = 0;

% amplitude of pattern
bi = [1,0,0];
% normalize
bi = bi/sqrt(sum(bi.^2));
b0 = sum(bi);
sb = 2;
% frequency (wave number) of pattern
fi = [1,2,3];

% length of domain for integration
L  = 1e3;
% scale of gaussian window were pattern is restricted to
Lw = 60;

% pattern function
mode = 'positive'
wfun = @(x) exp(-1/2*(x/Lw).^2);
switch (mode)
case {'positive'}
%	yfun = @(x) sqrt(2)*(1+cos(2*pi*x*ki))*bi'.*wfun(x); %normpdf(x/Lw);
	yfun = @(x) sqrt(2)*sb^2*(b0+cos(2*pi*x*fi)*bi').*normpdf(x/Lw);
otherwise
	yfun = @(x) a*sqrt(2)*cos(2*pi*x./ll).*normpdf(x/Lw);
	% note : the strictly positive pattern concentrates half its energy in the mean
end
%yfun = @(x)a/(ll)*(1+cos(2*pi*x/ll)).*(x>-1/2*Lw).*(x<=1/2*Lw);

% transect coordinates
x  = innerspace(-L/2,L/2,n)';
dx = x(2)-x(1);

% samples pattner
yr = se.*randn(n,1);
yr = yr-min(yr);
yr = yr.*wfun(x);
y = yfun(x) + yr;
%y = y-mean(y);

% root of spectral energy
% note that the spectral density for the continous transform is normalized
% by the spectral energy
% note that the spectral energy diverges when the pattern is not restricted
% to a square-integrable window
rse = sqrt(sum(y.^2)*dx)

% restriction to a window spreads the spectral energy over a frequency range,
% so evaluate around the central frequency
fc = mean(fi);
fx = linspace(0,3*fc,nf);
df = fx(2)-fx(1);
% continuous fourier transform
f = cos(2*pi*x*fx)'*y*dx + 1i*sin(2*pi*x*fx)'*y*dx;
% spectral density
S = 2/(rse^2)*(f.*conj(f));
% forward transform for R
x_ = x(x>=0);
R = cos(2*pi*x_*fx)*S*df + 1i*sin(2*pi*x_*fx)*S*df;

% integrate spectral energy, check that this values is 1
IS = sum(S.*df)

%bmu = mean(y);

figure(1)
clf
subplot(2,2,1)
plot(x,y)
hline([0.5*max(y)])
subplot(2,2,2)
cla();
hold on
plot(fx,S,'.-');
hold on

% compare to analytic spectrum
S_ = spectral_density_wperiodic(fx,Lw,fi,bi);
plot(fx,S_);
IS(2) = sum(S_*df)

fx_per = fourier_axis(x);
S_per = periodogram(y,x(end)-x(1),length(y));
plot(fx_per,S_per,'-');
IS(3) = sum(S_per)*(fx_per(2)-fx_per(1))
legend('quad','anal','per')

%S_ = 1/(sqrt(pi)*Lw*(sum(bi.^2)+2*sum(bi).^2))*(f.*conj(f));
%rse_ = sqrt(2*Lw*sqrt(pi));
%rse_ = 1*sqrt(sum(S_.*dk));
%S_ = 2/rse_.^2*S_;

vline(fi);
ylim([0,1.1*max(real(S))]);
xlim([0,1+max(fi)])

R(:,2) = autocorrelation_wperiodic(x_,x_(1),Lw,fi,bi);
%R(:,2) = wfun((x_-x_(1))/sqrt(2)).*(2*b0^2+cos(2*pi*(x_-x_(1))*ki)*(bi).^2')/(1+2*b0^2);

subplot(2,2,3);
plot(x_,[real(R)]); %,imag(R)]);

R_ = autocorr_fft(y);
hold on
plot(x-x(1),R_);


figure(2)
clf
p = 0.25;
o = 0.0;
pw = 0.0;
msk = (x>-p*Lw & x<(p*Lw+o));
w = tukeywin(sum(msk),pw);

%subplot(2,2,4)
x_ = x(x>= 0 & x < 2*p*Lw);
y_ = y(msk);
y_ = y_ - mean(y_);
y_ = w.*y_;
% continuous fourier transform
f = cos(2*pi*x(msk)*fx)'*y_*dx + 1i*sin(2*pi*x(msk)*fx)'*y_*dx;
rse = sqrt(sum(y_.^2)*dx);
% spectral density
S = 2/(rse^2)*(f.*conj(f));
R = cos(2*pi*x_*fx)*S*df + 1i*sin(2*pi*x_*fx)*S*df;
l = length(R);
R = l./(l-(1:l)'+1).*R/R(1);

subplot(2,2,1)
plot(x(msk),y_)

subplot(2,2,2)
plot(fx,S)
hold on
fx_ = fourier_axis(x(msk));
plot(fx_,periodogram(y_,2*p*Lw,length(y_)));
xlim([0,1+max(fi)]);

subplot(2,2,3)
%hold on
plot(x_,[real(R)]); %,imag(R)]);
R_ = autocorr_fft(y_);
hold on
x_ = x(msk);
plot(x_-x_(1),R_);
ylim([-1.5,1.5])


% x=linspace(-5,5,1e3)'; w=3; s=1/(2*pi*w); y=[normpdf(x,0,s).^2, 1/(s*2*sqrt(pi))*normpdf(x,0,s*1/sqrt(2))], sum(y)*(x(2)-x(1))
% x=linspace(-5,5,1e3)'; plot([normpdf(x).^2, 1/(2*sqrt(pi))*normpdf(x,0,1/sqrt(2))])


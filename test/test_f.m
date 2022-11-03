% 2021-09-14 14:56:26.733152500 +0200

 k = 10*[1,2,3];
 k = 10;
 bi = [1,3,5];
 bi = 1;
 bi = bi./sqrt(sum(bi.^2));
 se = 10;
 sb = 3;
m = 1e4;
n = 1e3;
L = 1;

% y=sqrt(2)*a*cos(20*pi*x);
% f=2*fft(y)/n;
% S=0.5/var(y)*abs(f).^2;
% fx=fourier_axis(x);
% df=fx(2)-fx(1);
% S_ = periodogram(y,1);
% max(S), sum(S(end/2:end))*df, std(y)

x = L*(0:n-1)'/n;
b = sqrt(2)*sb*cos(2*pi*x*k)*bi';

% random perturbation
b = b + se*randn(n,m);
b = b - mean(b);

f = fft(b)/n;
S = f.*conj(f);

df = 1/L;
before = false;%true;
if (before)
scale = 1./(sum(S)*df);
S = scale.*S;
end


S_mu = mean(S,2);
S_sd = std(S,[],2);

if (~before)
scale = 1/(sum(S_mu)*df);
S_mu = scale*S_mu;
S_sd = scale*S_sd;
end
fx = fourier_axis(x);
clf
plot(fx,[S_mu,S_sd],'o');
hold on
den = (1/4*sb^2 + 1/2*se.^2)*df
S_k = 1/den*(1/4*sb.^2*bi.^2+1/2*se^2)
plot(k,S_k,'*')
hline(1/n*se^2/(sb^2 + se^2))
plot(k,sqrt(1/16*1/n*2*se*bi),'*');


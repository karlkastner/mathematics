% Tue 18 Jan 11:43:48 CET 2022

L = 1e3;
n = 1e5;

fc = 1;
p = [1,2,3];

IS = spectral_density_lorentzian_scale(fc,p)

x = linspace(0,L,n)';
fx=fourier_axis(x);

S = spectral_density_lorentzian(fx,fc,p,-1);

%df  = (fx(2)-fx(1))
df = 1/L;
fdx = fx>=0;
IS_ = nansum(S(fdx,:),1)*df




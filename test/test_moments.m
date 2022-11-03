% only mean and std match, higher moments not!
n=1e6;
 L = 1e3;
 x = linspace(0,L,n);
 f=fourier_axis(x);
 fdx=f>0;
 S = exp(-abs(f));
 wmean(S(fdx),f(fdx)), wstd(S(fdx),f(fdx)), plot(f,S), Sh = (1/2*chi2rnd(2,n,1)).*S;
 wmean(Sh(fdx),f(fdx)), wstd(Sh(fdx),f(fdx));
 wskew(Sh(fdx),f(fdx)), wkurt(Sh(fdx),f(fdx));
 plot(f,[S,Sh])
kurtosis(1/2*chi2rnd(2,n,1))

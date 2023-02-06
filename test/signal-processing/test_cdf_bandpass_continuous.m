fx = linspace(0,1e1,1e2)';
fc = 2;
p  = 3;
S  = spectral_density_bandpass_continuous(fx,fc,p);
cdf = cumsum(S).*(fx(2)-fx(1));
cdf(:,2) = cdf_bandpass_continuous(fx,fc,p);
%cumsum(S).*(fx(2)-fx(1));

plot(fx,cdf)


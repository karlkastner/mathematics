fx = linspace(0,1e1,1e2)';
fc = 2;
sx  = 0.3;
S  = spectral_density_brownian_phase(fx,fc,sx);
cdf = cumsum(S).*(fx(2)-fx(1));
cdf(:,2) = cdf_brownian_phase(fx,fc,sx);
%cumsum(S).*(fx(2)-fx(1));

plot(fx,cdf)


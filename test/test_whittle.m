n = 100;
L = 1;
df = 1/L;
x = linspace(0,L,n)';
fx = fourier_axis(x);
p = 1; 
S=spectral_density_bp(fx,[1:21],L,n,p,'f');
 S=2*S./(sum(S)*df);
 plot(S);


[par,S_] = fit_spectral_density(fx,S(:,11),[11,p]/2,L,'f','bp') 
[par,S__] = fit_spectral_density(fx,S(:,11),[11,p],L,'f','lorentzian') 

clf
subplot(2,2,1)
plot([S(:,11),S_,S__])

 S=S(2:end,:);
 LL = sum(log(S)) + sum(S(:,(end+1)/2)./S);
subplot(2,2,2)
 plot(1:21,LL);
 ylim(quantile(LL,[0,0.5]));
 vline(11);



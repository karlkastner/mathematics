% 2022-04-26 15:42:01.216084769 +0200

m  = 1000;
m  = 100;
dx = 0.02;
L  = 100; 

f0 = 1;
s = 0.3;

% integral
I = quad(@(fx) spectral_density_brownian_phase(fx,f0,s),0,10000)

[xi,yi,phi] = fourier_random_phase_walk(1,0,f0,s,L*m,dx);
fx=fourier_axis(xi);
yi = yi-mean(yi);
S = periodogram_bartlett(yi,L*m,m,length(xi),[],false);

S(:,2) = spectral_density_brownian_phase(fx,f0,s);
[fc,Sc] = spectral_density_brownian_phase_mode(f0,s)         
%S(:,3) = spectral_density_brownian_phase(fx,f0,s*1.1);
%S(:,4) = spectral_density_brownian_phase(fx,f0,s/1.1);

fdx = fx>0;
clf
plot(fx(fdx),S(fdx,:))
%vline(f0)
hold on
plot(fc,Sc,'o')
vline(fc,'color','r')
axis auto

% Sat 11 Jun 15:20:04 CEST 2022
sy=0.4;
n=1e6;
fy = 2*linspace(-1,1,n)';
S=spectral_density_brownian_phase_across(fy,sy);
plot(fy,S)
s=sum(S*(fy(2)-fy(1)))

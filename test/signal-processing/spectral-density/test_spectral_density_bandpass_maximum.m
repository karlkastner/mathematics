fc = 2;
order = 3;
f = linspace(0,10*fc,1e3);

[S,IS] = spectral_density_bandpass_continuous(f,fc,order);
IS
[Sc,IS] = sd_bandpass_continuous_max(fc,order);
IS
clf
plot(f,S)
hold on
plot(fc,Sc,'o')

p = sd_bandpass_max2par(fc, Sc)
 

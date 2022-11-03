% Tue 18 Jan 11:08:00 CET 2022

f0 = 2;
%p = 1./[2,4,9];
%p = [2,4,9];
p = [1,2,4];
%fc = matlabFunction(fc);
%Sc = matlabFunction(Sc);
L = 1e2;
n = 1e5;

fc = f0;
[Sc] = spectral_density_lorentzian_max(f0,p);


fx = linspace(0,L,n)';
clf();
[S,IS] = spectral_density_lorentzian(fx,f0,p); %.^(-1/2)))

plot(fx,S);
hold on
plot(fc,Sc,'*')
vline(fc)

%fc = fc.';
%Sc_ = spectral_density_lorentzian(fc*[0.999 1.001],f0,p',false);
%Sc' ./ Sc_ < 1



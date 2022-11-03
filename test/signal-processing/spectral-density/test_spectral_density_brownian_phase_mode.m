% 2022-01-06 17:20:01.004124047 +0100
% note : modes S^p are not defined for p>0

function fail = test_spectral_density_brownian_phase_mode()

syms f f0 s
[fc,Sc] = spectral_density_brownian_phase_mode(f0,s)
S = spectral_density_brownian_phase(fc,f0,s)

% value of Sc
res1fun = matlabFunction(S-Sc);
f0_ = 1.1;
s_ = 0.3;
res(1) = rms(res1fun(f0_,s_));

% check maximum, dS/df = 0
S = spectral_density_brownian_phase(f,f0,s);
dS_df  = diff(S,f);
dS_df = simplify(subs(dS_df,f,fc));

resfun2 = matlabFunction(dS_df,'Vars',{'f0','s'});
res(2) = rms(resfun2(f0_,s_));

fail = rms(res) > 1e-4;

if (0)
f0 = 2;
s = [1.4,4,9];

L = 1e2;
n = 1e5;

fx = linspace(0,L,n)';
clf;
[S] = spectral_density_brownian_phase(fx,f0,p); %.^(-1/2)))
plot(fx,S); %(fx,f0,p.^(1)))
hold on
plot(fc,Sc,'*')
vline(fc)

fc = fc.';

Sc_ = spectral_density_brownian_phase(fc*[0.999 1.001],f0,p',false);
Sc' ./ Sc_ < 1
end

end

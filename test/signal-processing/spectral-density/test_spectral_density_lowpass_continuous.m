 syms f fc p;
 [S,S1] = spectral_density_lowpass_continuous(f,fc,p,true), I=int(S1,f), subs(subs(I,fc,2),f,inf)

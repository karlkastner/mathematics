% Mon 26 Jul 13:22:12 CEST 2021
syms f r0 L n p
S = spectral_density_bp(f*n,r0*L,L,n,1/2,false)
solve(S-1,r0)

syms f r0 L n p S0;
S=spectral_density_bp(f*n,r0*L,L,n,p,false);
solve(S-S0,p)


% standard triangle

syms chi psi

phi_1 = 1 - chi - psi
phi_2 = chi - psi

subs(subs(int(phi_1^2,chi),chi,psi),psi,1)
subs(subs(int(phi_1*phi_2,chi),chi,psi),psi,1)



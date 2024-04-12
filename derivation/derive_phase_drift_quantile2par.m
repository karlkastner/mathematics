% 2024-01-27 16:17:49.124952744 +0100

syms f0 s s2 a2 f1 f2 a1 C1 C2; pi_=sym(pi)
f = [f1,f2]
C=phase_drift_cdf(f,f0,s); 

% C = atan2(2*f*f0*s^2*pi, f0^2 - f.^2 + f0^2*s^4*pi^2)/pi

eq1 = (a1*(f0^2 - f1^2 + f0^2*s2^2*pi_^2) - 2*f1*f0*s2*pi_)
eq2 = (a2*(f0^2 - f2^2 + f0^2*s2^2*pi_^2) - 2*f2*f0*s2*pi_)

f01 = solve(eq1,f0)
f02 = solve(eq2,f0)

f01_ = (f1*(a1^2*pi*s2^2 + a1^2 + pi^2*s2^2)^(1/2) + pi*f1*s2)/(a1*pi*s2^2 + a1)
f02_ = (f2*(a2^2*pi*s2^2 + a2^2 + pi^2*s2^2)^(1/2) + pi*f2*s2)/(a2*pi*s2^2 + a2)

s2_ = solve(f01==f02,s2)

s21 = solve(eq1,s2)
s22 = solve(eq2,s2)

solve(s21(1)==s22(1),f0)


% Tue 30 Jan 14:19:51 CET 2024

f0 = 1.1
s  = 0.3

f = [0.5*f0,3*f0];
C = phase_drift_cdf(f,f0,s)

s2 = s*s;
f1 = f(1);
f2 = f(2);
a1 = atan(pi*C(1))
a2 = atan(pi*C(2))

f0_ = [
 (f1*(a1^2*pi*s2^2 + a1^2 + pi^2*s2^2)^(1/2) + pi*f1*s2)/(a1*pi*s2^2 + a1)
-(f1*(a1^2*pi*s2^2 + a1^2 + pi^2*s2^2)^(1/2) - pi*f1*s2)/(a1*pi*s2^2 + a1)
]
f0 - f0_

[f0_,s_] = phase_drift_quantile2par(f(1),f(2),C(1),C(2))
f0_ - f0
s_ - s

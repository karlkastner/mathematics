% Thu 20 Jul 14:45:20 CEST 2017

syms r1 r2 a1 a2
a = [-r1/(r2 - 1), -(r1^2 - r2^2 + r2)/(r2 - 1)], % term one and 2 of acf (term 0 is 1)


r2_ = solve(a1-a(1),r2)
r1_ = solve(a1-a(1),r1)

r1 = solve(subs(a2-a(2),r2,r2_),r1)
r2 = solve(subs(a2-a(2),'r1',r1_),r2)



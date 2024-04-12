% 2024-01-15 11:10:10.171248642 +0100

syms CC_ C_ a b mu_ s2_
mu = a/(a+b)
s2 = a*b/((a+b)^2*(a+b+1))
C = b/((a+b)*(a+b+1))
CC = b/(a*(a+b+1))

solve(C-C_,b)

%s2/mu2 = b./(a*(a+b+1))
%s2_/mu2_*(a+ = b./(a*(a+b+1))
%solve(mu-mu_,s2_-s2_,[a,b])

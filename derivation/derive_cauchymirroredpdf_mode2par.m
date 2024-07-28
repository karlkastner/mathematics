syms f fc s f0 Sc;
syms s2 positive
S = cauchymirroredpdf(f,f0,s);

if (0)
[n,d]=numden(S);
 S = n/d;
 [n,d] = numden(diff(S,f));
 fcs = solve(n,f), fcfun = matlabFunction(fcs);
 fcfun(-1,1), solve(fcs(2)-fc,f0), solve(fcs(3)-fc,f0)
end
eq = subs(S,f,fc) - subs(S,f,-fc)
sol_f0p = solve(subs(S,f,fc) - Sc,f0)
sol_f0n = solve(subs(S,f,-fc) - Sc,f0)

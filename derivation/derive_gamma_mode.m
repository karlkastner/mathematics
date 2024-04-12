syms a b x; x0=(a-1)/b;  f = b^a/gamma(a)*x^(a-1)*exp(-b*x), eq=simplify(diff(f,x)./f), solve(eq,x) 

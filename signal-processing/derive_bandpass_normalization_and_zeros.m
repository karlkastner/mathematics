% 2021-07-12 14:53:47.505426227 +0200

syms f f0 f0_ fmax;
%F=int(1./(1-2*f0*cos(pi*f)+f0^2),f), subs(F,f,1), subs(F,f,0) 

if (0)
 h = (1-2*(f0/fmax) + (f0/fmax)^2)./(1-2*f0/fmax*cos(pi*f/fmax)+(f0/fmax)^2);
 h = h*(1-h);
 dh = diff(h,f);
 s = solve(dh,f)
end

%simplify(solve(f0_ - s(2),f0))
 syms f r0 f0 L n;
 s = spectral_density_lp(f,r0,L,n,1/2,'rho'),
 s = solve(diff(s*(1-s),f),f)
 simplify(solve(s(2) - f0,r0))
 simplify(solve(s(3) - f0,r0))


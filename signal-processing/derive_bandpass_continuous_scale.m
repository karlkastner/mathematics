% 2021-07-12 14:53:47.505426227 +0200

if (0)
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

end

% wolfram : integrate 1/2pi*( (2*k)/(1 + k^2) )^(2*p) dk
% = 1/2pi*4^p k^(2p + 1) / (2p + 1) 2F1(2p, p+1/2,p+3/2,-k^2) 

kc = 1/3;

S = @(k,p) ((2*k/kc)./(1 + (k/kc).^2)).^(2*p)

L = 1e5;

Ifun = @(k,p) kc/(2*pi) * 4^p.*k.^(2*p+1) ./ (2*p+1).*hypergeom([2*p,p+1/2],p+3/2,-k.^2)

Ifun2     = @(k,p) kc/(2*pi) * 4^p.*k.^(2*p+1) ./ (2*p+1).*(1 + k.^2).^(-p-1/2).*hypergeom([p+3/2-2*p,p+1/2],p+3/2,-k.^2./(-k.^2-1))
Ifun_inf  = @(p) kc/(2*pi) * 4^p./(2*p+1).*hypergeom([p+3/2-2*p,p+1/2],p+3/2,1)
Ifun_inf  = @(p) kc/(2*pi) * 4^p./(2*p+1).*hypergeom([3/2-p,p+1/2],p+3/2,1)
Ifun_inf2 = @(p) kc/(2*pi) * 4^p./(2*p+1).*gamma(p+3/2).*gamma(p-1/2)./gamma(2*p)./gamma(1);
Ifun_inf2 = @(p) kc/(2*pi) * 4^p./(2*p+1).*gamma(p+3/2).*gamma(p-1/2)./gamma(2*p);
Ifun_inf2 = @(p) kc/(2*pi) * 4^p./(2*p+1).*(p+1/2).*gamma(p+1/2).*gamma(p-1/2)./gamma(2*p);
Ifun_inf2 = @(p) kc/(2*pi) * 2.^(2*p-1).*gamma(p+1/2).*gamma(p-1/2)./gamma(2*p);
Ifun_inf2 = @(p) kc/(2*pi) * 2.^(2*p-1)./(p-1/2)*gamma(p+1/2)^2./gamma(2*p);
% legendre duplication
Ifun_inf2 = @(p) kc/(2*pi) * 2.^(2*p-1)./(p-1/2)*(2^(1-2*p)*sqrt(pi)*gamma(2*p)./gamma(p))^2./gamma(2*p);
Ifun_inf2 = @(p) kc/(2*pi) * 2.^(1-2*p)./(p-1/2)*pi*(gamma(2*p)./gamma(p)^2);
Ifun_inf2 = @(p) kc/(2*pi) * 2.^(2-2*p)./(2*p-1)*pi*(gamma(2*p)./gamma(p)^2);
Ifun_inf2 = @(p) kc * 2.^(1-2*p)./(2*p-1)*(gamma(2*p)./gamma(p)^2);
Ifun_inf2 = @(p) kc *1/sqrt(pi)*1./(2*(p-1/2))*(gamma(p+1/2)./gamma(p));
Ifun_inf2 = @(p) kc *1/(2*sqrt(pi))*gamma(p-1/2)./gamma(p);
%Ifun_inf2 = @(p) kc/(2*pi) * 4.^(1-p)./(2*p-1)*pi*p*binom(2*p-1,p-1);

pp = innerspace(0,4)';

I = [];
TOL = 1e-14;
for idx=1:length(pp)
	p = pp(idx);
	I(idx,1) = 1/(2*pi)*quadl(@(k) S(k,p),0,L,TOL);
	try
	I(idx,2) = Ifun(L,p);
%	I(idx,3) = Ifun2(L,p);
	I(idx,3) = Ifun_inf(p);
	I(idx,4) = Ifun_inf2(p);
	catch
		e
	end
end
subplot(2,2,1)
plot(pp,I)
ylim([-1,20])
%semilogy(pp,[real(I),imag(I)])
xlabel('p');
ylabel('I');
subplot(2,2,1)

syms k p
S = ((2*k/kc)./(1 + (k/kc).^2)).^(2*p);
pp = [0.5:0.5:3];
I = [];
for idx=1:length(pp)
	I(idx,1) = 1/(2*pi)*double(int(subs(S,p,pp(idx)),k,0,inf))
	I(idx,2) = Ifun_inf2(pp(idx))
end


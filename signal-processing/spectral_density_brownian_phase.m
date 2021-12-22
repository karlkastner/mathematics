% Wed  1 Dec 18:55:18 CET 2021
function S = spectral_density_brownian_phase(fx,fc,s)
	%p  = pi*s.^2;
	p = pi./s.^2;
	%S = (p*(k.^2 + p.^2 + kc.^2))/( (k.^2 + p.^2 - kc.^2).^2  + 4*kc.^2*p.^2)
	S = (fx.^2./fc.^2 + p.^2 + 1)./( (fx.^2./fc.^2 + p.^2 - 1).^2  + 4*p.^2);
	if (isnumeric(fx))
	if (0)
		fmax = max(fx);	
		S = p./(f0*(atan(p*(fmax/f0 - 1)) + atan(p))).*S;
	else
		I = spectral_density_area(fx,S);
		S = S./I;
	end
	else
		S=S./(fc*s^2);
	end
end % spectral_density_brownian_phase

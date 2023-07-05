% Sat  2 Jul 15:29:09 CEST 2022
% note : this still excludes the small side lobe near zero,
%        and is thus only valid for small s
function [p,pa] = brownian_phase_patch_size_distribution(l,fc,s)
	l = cvec(l);
	omega = 2*pi*fc;
	%c1 = normcdf(0,omega*l,s*sqrt(l)) + (1-normcdf(pi,omega*l,s*sqrt(l)));
	%p1 = normcdf(-omega*l./(s*sqrt(l))) + (1-normcdf((pi-omega*l)./(s*sqrt(l))));
	%p2 = 1 - normcdf((omega*l-pi)./(s*sqrt(l))) - exp(2*pi*omega)*(1-normcdf((omega*l+pi)./(s*sqrt(l))));
	p1 = normcdf(-omega*l./(s*sqrt(l))) + (normcdf((omega*l-pi)./(s*sqrt(l))));
	p2 = 1 - normcdf((omega*l-pi)./(s*sqrt(l))) - exp(2*pi*omega)*(normcdf((-omega*l-pi)./(s*sqrt(l))));
	p  = p1.*p2;
	% approximation
	pa = (normcdf((pi-omega*l)./(s*sqrt(l)))).*(normcdf((omega*l-pi)./(s*sqrt(l))));
	p(l==0)  = 0;
	pa(l==0) = 0;
	dl = diff(l);
	p  = p./(sum(mid(p).*dl));  
	pa = pa./(sum(mid(pa).*dl));  
end


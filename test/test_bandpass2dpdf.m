%fun = @(fr,fc,p) 
syms fr fc p x positive 
fr = x
%p = 4;
%S = (1./(1 + (fr/fc).^2).^2.*( 1 - 1./(1 + (fr/fc).^2).^2)).^p;
%S = (1./(1 + (fr/fc).^2).*( 1 - 1./(1 + (fr/fc).^2))).^(2*p);

xSfun = @(fr,fc,p) 2*pi*fr.*(4*1./(1 + (fr/fc).^2).*( 1 - 1./(1 + (fr/fc).^2))).^(2*p);

%p=2;
% x=1e4;
% a=2*p;
% b=p+1;
% c=p+2;
% (-z)^(p)*(-z)^(-2*p+1)*hypergeom([c-a,c-b],c,z)/(2*p+2), xSfun = @(x,p_) x.*(x.^2./(1+x.^2).^2).^p_;

fc = 0.5;
p = 1.5;

integral(@(x) xSfun(x,fc,p),0,inf)


if (1)
	% normalization of Sxy bandpass 
	pi*fc^2*2^(4*p)/(2*p+1)*gamma(2*p+2)*gamma(2*p-1)/(gamma(4*p));
end

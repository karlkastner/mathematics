% Mon 20 Sep 10:21:10 CEST 2021
% spectral density of a periodic function restricted by a gaussian window
% S = |f|^2
% y = exp(-1/2*(x/Lw)^2)*(b0 + bi*sum(2*pi*fi*x))
function S = spectral_density_wperiodic(fx,Lw,fi,bi)
	bi = bi/sqrt(sum(bi.*bi));
	b0 = sum(bi);
	s_w = 1/(2*pi*Lw);
	f  = (        normpdf(cvec(fx),fi,s_w)*bi' ...
	       + 2*b0*normpdf(cvec(fx),0,s_w) );
	S  = 1/(sqrt(pi)*Lw*(1+2*b0.^2))*(f.*conj(f));
end


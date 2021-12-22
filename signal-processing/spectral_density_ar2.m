% Fri 24 Sep 19:27:17 CEST 2021
% Bartlett, Berkley Lecture notes
function S = spectral_density_ar2(fx,c)
	p = roots(c);
	fx = 0.5*fx/max(fx);
	S = 1./(abs(exp(-2i*pi*fx) - p(1)).*abs(exp(-2i*pi*fx)-p(2))).^2;
	df = fx(2)-fx(1);
	fdx = fx>= 0;
	S = 2*S/(sum(S(fdx))*df);
end


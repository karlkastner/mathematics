% Mon 20 Sep 10:17:01 CEST 2021
% autocorrelation of a periodic function restricted by a gaussian window
% R(x') = E[(y(x)-bary)(y(x+x')-bar y)]/sigma_y^2
% y = exp(-1/2*(x/Lw)^2)*(b0 + bi*sum(2*pi*fi*x))
function R = acf_periodic_windowed(x,x0,Lw,fi,bi)
	bi = bi/sqrt(sum(bi.*bi));
	b0 = sum(bi);
	x = x-x0;
	R = exp(-1/4*(x.*x/(Lw*Lw))) ...
	    .*(2*b0*b0 + cos(2*pi*x*fi)*(bi.*bi).')/(1+2*b0*b0);
end



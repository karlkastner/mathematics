% Sat 11 Aug 19:46:02 CEST 2018
% Karl Kastner, Berlin 
%
% continous fourier coefficients
%
% TODO what does this function actually do again?
function [a,b] = hyperbolic_fourier_box(n,L,w,lval,rval)
	a =  (lval - rval.*exp(-(n*w)/L))./(exp((n*w)/L) - exp(-(n*w)/L));
	b = -(lval - rval.*exp((n*w)/L))./(exp((n*w)/L) - exp(-(n*w)/L));
	if (~issym(n))
		a(n == 0) = 1;
		b(n == 0) = 0;
	end
end


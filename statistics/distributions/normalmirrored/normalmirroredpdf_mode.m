% Mon 13 Feb 13:42:27 CET 2023
function [fc,Sc] = normalmirrored_mode(fm,s)
	% d/df ( 
	%syms fm s fc
	%S = exp(-1/2*(f-fm)^2/s^2) + exp(-1/2*(f+fm)^2/s^2)
	%$S = exp(-1/2*(f^2 - 2*f*fm + fm^2)/s^2) + exp(-1/2*(f+fm)^2/s^2)
	%dS/df = (f-fm)*exp(-2*f*fm) + (f+fm)*exp(+2*f*fm) = 0
	%fun = @(f) (f-fm).*exp(-0.5*(f-fm).^2/s^2) + (f+fm).*exp(-0.5*(f+fm).^2/s^2);
	fun = @(f) (f-fm).*exp(-0.5*(f-fm).^2/s^2) + (f+fm).*exp(-0.5*(f+fm).^2/s^2);
	fc = fzero(fun,fm);
	Sc = normalmirroredpdf(fc,fm,s);
end

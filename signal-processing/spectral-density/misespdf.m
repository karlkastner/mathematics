% Sun 29 Jan 10:13:02 CET 2023
function f = misespdf(x,mu,k)
%	den = besseli(0,k);
	f = exp(k.*cos(x-mu))./(2*pi*besseli(0,k));
	f_large = exp(k.*cos(x-mu) - besseliln_large_x(0,k))./(2*pi); 
	fdx = ~isfinite(f);
	f(fdx) = f_large(fdx);
end



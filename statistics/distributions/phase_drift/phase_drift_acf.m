% Wed  1 Dec 10:49:20 CET 2021
function R = acf_brownian_phase(x,fc,s)
	kc = 2*pi*fc;
	R = cos(kc.*x).*exp(-pi*kc*s^2*abs(x));
end


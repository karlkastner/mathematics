% Wed  1 Dec 10:49:20 CET 2021
function R = phase_drift_acf(x,f0,s)
	k0 = 2*pi*f0;
	R = cos(k0.*x).*exp(-pi*k0*s^2*abs(x));
end


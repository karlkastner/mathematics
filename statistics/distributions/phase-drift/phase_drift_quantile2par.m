% Tue 30 Jan 14:48:37 CET 2024
function [f0,s] = phase_drift_quantile2par(f1,f2,C1,C2)
	a1 = tan(pi*C1);	
	a2 = tan(pi*C2);
	tmp = sqrt( -a1^2*a2^2*(f1^2 - f2^2)^2 + 4*(a1^2+a2^2)*f1^2*f2^2 - 4*a1*a2*f1*f2*(f1^2 + f2^2))
	s2 = (a1*a2*(f1^2 - f2^2))/(pi*tmp);
	s2 = abs(s2);
	s  = sqrt(s2);
	f0 = tmp/(2*a1*f2 - 2*a2*f1);
	f0 = abs(f0);
	%f0 =  (- a1^2*a2^2*f1^4 + 2*a1^2*a2^2*f1^2*f2^2 - a1^2*a2^2*f2^4 + 4*pi*a1^2*f1^2*f2^2 - 4*pi*a1*a2*f1^3*f2 - 4*pi*a1*a2*f1*f2^3 + 4*pi*a2^2*f1^2*f2^2)^(1/2)/(pi^(1/2)*(2*a1*f2 - 2*a2*f1));
end


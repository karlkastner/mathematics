% 2018-02-02 17:53:28.117627900 +0100
%
% cutoff frequency of the ar1 filter
%
% H = a / (1-(1-a)z^-1)
% a = 1-exp(-omega_c * T)
function [omega_c, T_c] = ar1_cutoff_frequency(rho,Ts)
	if (nargin()<2)
		Ts = 1;
	end
%	T_c = acos((rho^2/2 - (2^(1/2)*(rho^2 - 2*rho + 1))/4 + 1/2)/rho)/(2*pi);
% 	-acos((rho^2/2 - (2^(1/2)*(2*rho^2 - 4*rho + 2))/4 + 1/2)/rho)/(2*pi)
%	p = sqrt(0.5);
	p = 0.5;
	% -acos((p + 2*rho + p*rho^2 - rho^2 - 1)/(2*p*rho))/(2*pi)
	f_c = acos((p + 2*rho + p*rho^2 - rho^2 - 1)/(2*p*rho))/(2*pi);
	% fimplified
	if (1)
	f_c(2) = ((1 - p)^(1/2)*(1 - rho))/(2*pi*p^(1/2))
	% which for p=1/2 is approximately
	% f_c \approx (1-rho)/(2*pi) \approx -log(rho)/(2*pi)
	end

%	f_c =  acos((rho^2/2 - (p*(2*rho^2 - 4*rho + 2))/4 + 1/2)/rho)/(2*pi);
%	f_c = 1/
	T_c = 1./f_c;
	omega_c = 2*pi.*f_c; %./T_c;
	
if (0)
	% y = s(x + rho yold) = ax + (1-a)yold
	s = (1-rho);
	a = 1-s*rho
	omega_c = -log(1-a)/Ts;
	f_c = omega_c/(2*pi);
	T_c = 1/f_c;
end
end


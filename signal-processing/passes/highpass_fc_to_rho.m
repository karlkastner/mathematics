% Wed 12 Jan 09:49:48 CET 2022
% cutoff frequency of the lowpass filter
% S(fc)^p = 0.5
function rho = highpass_fc_to_rho(fc,p,dx)
	if (nargin()<2)
		p = 1;
	end
	a   = 0.5^(1/(2*p));
	f = fc;
	b = cos(2*pi*dx*f);
	rho = (1 + a*b - b + sin(pi*dx*f)*sqrt(2*(1 - a)*(1 + a - b + (a*b))))/a;

end


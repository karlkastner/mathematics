% Wed 12 Jan 09:49:48 CET 2022
% cutoff frequency of the lowpass filter
% S(fc) = 0.5
function rho = lowpass_fc_to_rho(fc,p,dx)
	if (nargin()<2)
		p = 1;
	end
	a   = 0.5^(1/(2*p));
	b   = a*cos(pi*dx*fc).^2;
	rho =  (1 + a - 2*b + 2*sqrt((b - 1).*(b - a)))./(1 - a);
end


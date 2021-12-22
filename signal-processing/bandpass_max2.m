% Thu  9 Sep 15:06:54 CEST 2021
% note - with the current implementation bp max is located at f0!
function fmax = bandpass_max2(f0,p,L,n)
% 	fmax(1) = -(n*(pi - acos(cos((2*pi*L*f0)/n) - 2^(1/(2*p))*cos((2*pi*L*f0)/n) + 2^(1/(2*p)) - 2) + 2*pi*k))/(2*L*pi);
%	fmax =  abs(n*(pi - acos(cos((2*pi*L*f0)/n) - 2^(1/(2*p))*cos((2*pi*L*f0)/n) + 2^(1/(2*p)) - 2) + 2*pi*k))/(2*L*pi);
%	fmax = abs(n*(pi - acos(cos((2*pi*L*f0)/n) - 2^(1/(2*p))*cos((2*pi*L*f0)/n) + 2^(1/(2*p)) - 2)))/(2*L*pi);
	a = 0.5.^(1./(2*p))
	fmax = abs(n*acos((cos((2*pi*L*f0)/n) - 1)/a - cos((2*pi*L*f0)/n) + 2))/(2*L*pi);

end


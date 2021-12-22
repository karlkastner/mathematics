% Sat 26 Jun 21:04:19 CEST 2021
function f0 = filter_r0_to_f0(r0,L,n)
	%fmax)
	if (nargin()<2)
		fmax = 1;
	end
	f0 = fminsearch(@(f0) filter_f0_to_rho(f0,L,n)-r0,1);
 	%f0 = fmax*(((cos((pi*r0)./fmax) - 1).*(cos((pi*r0)./fmax) - 3)).^(1/2) - cos((pi*r0)./fmax) + 2);
end



% Sat 26 Jun 21:04:19 CEST 2021
function f0 = filter_rho_to_f0(rho0,dx)
	%fmax)
	if (nargin()<2)
		fmax = 1;
	end
	opt.Display = 'off';
	f0 = lsqnonlin(@(f0) filter_f0_to_rho(f0,dx)-rho0,1,[],[],opt);
 	%f0 = fmax*(((cos((pi*r0)./fmax) - 1).*(cos((pi*r0)./fmax) - 3)).^(1/2) - cos((pi*r0)./fmax) + 2);
end



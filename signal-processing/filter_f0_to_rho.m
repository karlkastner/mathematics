% Fri 16 Jul 17:12:41 CEST 2021
% function rho0 = filter_f0_to_rho(f0,L,n)
function rho0 = filter_f0_to_rho(f0,dx)
	if (nargin()<2)
		fmax = 1;
	end
%	omega0 = 2*pi*f0;
%	rho0 = 2 - cos(omega0*L/n) + (cos(omega0*L/n)^2 - 4*cos(omega0*L/n) + 3)^(1/2); 
 	rho0 = 2 - (cos(2*pi*dx*f0).^2 - 4*cos(2*pi*dx*f0) + 3).^(1/2) - cos(2*pi*dx*f0);
	if (~issym(rho0))
		fdx = (rho0 > 1);
		rho0_ =  2 + (cos(2*pi*dx*f0).^2 - 4*cos(2*pi*dx*f0) + 3).^(1/2) - cos(2*pi*dx*f0);
		rho0(fdx) = rho0_(fdx);
	end
end


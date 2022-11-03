% Fri 16 Jul 17:12:41 CEST 2021
% Karl Kästner, Berlin
%
%% correlation coefficient for the pth-order symmetric bandpass filter with
%% maximum at f0 (when rho_lp = rho_hp)
%
% note that this is independent of the filter order p
%
% function rho0 = bandpass_f0_to_rho(f0,L,n)
function rho0 = bandpass_f0_to_rho(f0,dx)
	if (nargin()<2)
		fmax = 1;
	end
	a    = cos(2*pi*dx*f0);
 	rho0 = 2 - a - sqrt(a.^2 - 4*a + 3);
	if (~issym(rho0))
		fdx = (rho0 > 1);
		rho0(fdx) = (4 - 2*a(fdx)) - rho0(fdx);
		%rho0(fdx) =  2 - a(fdx) + sqrt(a(fdx).^2 - 4*a(fdx) + 3);
		%rho0(fdx) = rho0_(fdx);
	end
end


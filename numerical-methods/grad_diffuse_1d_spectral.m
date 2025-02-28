% Wed  6 Nov 11:03:29 CET 2024
% compute dy(t+dt)/de, where y(t+dt) is determined (exaclty) by spectral methods:
% y(t+dt) = F^-1(exp(-e*dt*D2)*(F*y*F'))F'
% F^-1 exp(-2*pi*(kx+ky)*d*dt) F y F' F^-T                        
% F^-1 (-2 \pi (kx+ky)*dt) exp(-2*pi*(kx+ky)*d*dt*e) F y          
% F^-1 (-2 \pi (kx+ky)*dt) F y_ F' F^-T
function dz_de = grad_diffuse_1d_spectral(z,dt,e,n,L)

	[fx] = fourier_axis(L,n);
	ox = 2*pi*fx;
	fz = fft(z);
	fz = exp(   (-dt*e)*ox.*ox ...
		).*fz;
	s = (-(ox).^2*dt);
	%s(1)=0;
	dfz_de = s.*fz;
	dz_de  = ifft(dfz_de);
	dz_de  = dz_de;
end


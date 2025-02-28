% Wed  6 Nov 11:03:29 CET 2024
% compute dy(t+dt)/de, where y(t+dt) is determined (exaclty) by spectral methods:
% y(t+dt) = F^-1(exp(-e*dt*D2)*(F*y*F'))F'
% F^-1 exp(-2*pi*(kx+ky)*d*dt) F y F' F^-T                        
% F^-1 (-2 \pi (kx+ky)*dt) exp(-2*pi*(kx+ky)*d*dt*e) F y          
% F^-1 (-2 \pi (kx+ky)*dt) F y_ F' F^-T
function dz_de = grad_diffuse_2d_spectral(z,dt,e,n,L)
	znext = reshape(z,n);

	[fx, fy] = fourier_axis_2d(L,n);
	fy = fy';
	fz = fft2(znext);
	fz = exp(   ((-4*pi*pi*dt*e(1))*fx.*fx ...
		   + (-4*pi*pi*dt*e(1))*fy.*fy) ...
		).*fz;
	s = ((-4*pi*pi*dt)*(cvec(fx.*fx)+rvec(fy.*fy)));
	%s = 1./((-4*pi*pi*dt)*(cvec(fx.*fx))).*1./((-4*pi*pi*dt)*rvec(fy.*fy));
	%s(1,:) = 0;
	%s(:,1) = 0;
	dfz_de = s.*fz;
	dz_de  = ifft2(dfz_de);
	dz_de  = flat(dz_de);
end


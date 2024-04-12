function [dz] = derivative_2_2d_spectral(z,dx)
		dt = 1;
		a = [0,0];
		e = [1,1];
		n = size(z);
		L = n.*dx;
		[fx, fy] = fourier_axis_2d(L,n);
		fy = fy';
		fz = fft2(z);
		D2 =((-4*pi*pi*dt*e(1))*fx + 2i*pi*dt*a(1)).*fx ...
		  + ((-4*pi*pi*dt*e(2))*fy + 2i*pi*dt*a(2)).*fy ;
		%D2 = D2/(4*pi);
		%n
		%D2 = D2/(n(1)*n(2));
		%D2 = NaN.*D2;
		dfz = D2.*fz;
		dz = real(ifft2(dfz));
end
	

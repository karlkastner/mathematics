% 2022-05-01 11:18:18.894047311 +0200

L = 500;
n = 500;

x = linspace(0,L,n);
fx = fourier_axis(x);	
fr = hypot(fx,fx');



	par = [1/50,2]
	S = spectral_density_bandpass_2d(fr,par(1),par(2));
	
	y = randn(n,n);
	y = ifft2(sqrt(S).*fft2(y));
	
	Shat = abs(fft2(y)).^2;

	[par,Sp,resn] = fit_spectral_density_2d(fx,Shat);
	par
	Shat = trifilt1(Shat,3);
	[par,Sp2,resn] = fit_spectral_density_2d(fx,Shat);

	par

figure(1)
clf
plot([S(1,:)',Sp(1,:)',Sp2(1,:)',Shat(1,:)'])

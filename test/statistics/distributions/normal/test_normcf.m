% 2023-12-07 08:39:05.203651690 +0100 
 L=100;
	n=L^2;
	 fx = fourier_axis(L,n);
	 S = normpdf(fx,2,1);
	R=real(fft(S));
	R=R/R(1);
	 x=fourier_axis(fx);
	 R(:,2) = normcf(x,2,1);
	 plot(real(R))

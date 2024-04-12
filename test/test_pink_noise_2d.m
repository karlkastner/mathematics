% 2024-01-05 10:57:41.640394431 +0100
	fc = 0.1;
	figure(1);
	s = 0.5;
	clf();
	for idx=1 %[0.5,1,2]

	 n = 1000*[1,1];
	 L = idx*n;
	% [fx,fy,frr] = fourier_axis_2d(L*[1,1],n*[1,1]);
	 %S = 1./(frr+fc).^2;
	 [S,T,fx,fy,frr] = pink_noise_2d_pdf(n,L);
	x = fourier_axis(fx);
	R = ifft2(S);
%	 [R,x,y,rr] = pink_noise_2d_acf(n,L);
%	 S = 1./(1+(frr./fc).^2);
	 if (fc == 0)
	 S(1) = 0;
	 end

	 R_=ifft2(S);
	 R_=R_/R_(1);
%	 x=fourier_axis(fx);
%	 y=fourier_axis(fy);
%	 r = sqrt(x.^2+y'.^2);
	 %R2=pink2d_acf(r);


	 subplot(2,2,1)
	 plot(fftshift(fx),fftshift(S(1,:)));
	 hold on
	 subplot(2,2,2)
	 loglog(fftshift(fx),fftshift(S(1,:)));
	 hold on
	 subplot(2,2,3)
	 plot(x,[R(:,1),R_(:,1)]);
	 hold on
	 %hold on;
	 %plot(R2(1,:));
	 %x(2)-x(1), fx(2)-fx(1)
	end
	figure(2)
	[e] = pink_noise_2d(n,L);
	clf
	s_ = std(e,[],'all');
	e = s/s_*e;
	subplot(2,2,1)
	imagesc(1+e);
	axis equal
	axis tight
	subplot(2,2,2)
	e_ = exp(e);
	imagesc(e_);
	axis equal
	axis tight


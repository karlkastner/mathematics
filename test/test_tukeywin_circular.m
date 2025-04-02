 n = 1024;
 L = n;
 theta = 0.25*L;
 x = fourier_axis(1,n);
fx = fourier_axis(L,n);
 z = exp(-abs(hypot(x,x'))/theta);
 z_=(z-min(z,[],'all'));

zend = z(1,end/2);

wz=(z-zend).*tukeywin_circular(n*[1,1],0.25)+zend;

S=abs(fft2(z)).^2;
wS=abs(fft2(wz)).^2;

figure(1)
subplot(2,2,1)
plot([z(:,1),wz(:,1)])
subplot(2,2,2)
%plot(log([S(:,1),wS(:,1)]));
loglog(fx,log([S(:,1),wS(:,1)]));
 %z(1,:));
 %log(S(1,:))) %imagesc(log(S)) 


% Tue 24 Dec 15:50:39 +08 2019
% low and high pass of fourier transform:
% amplitude, frequency and phase remain the same
% -> multiplication with one large sine/cosine in real space
% low pass : tapered at ends (1 in centre, 0 at boundary)
% high pass: (1 at ends, 0 at center)

% -> does not drop exponentially, even with tukey window -> fit laplacian exp(-abs(dx))

L=2; x_ = linspace(0,L,1024); y=2*sin(2*pi*x_*10)'; figure(1); clf; plot(x_,y); fa=(abs(fft(y)))/length(x); max(fa), %plot(abs(fa)); max(fa)
plot([y,real(ifft(conv(fft(y),-[1,-2,1]/4,'same'))),real(ifft(conv(fft(y),[1,2,1]/4,'same')))])

L=100; x_ = linspace(0,L,5000); y=sin(2*pi*x_/L*500)'; figure(1); clf; plot(x_,y); fa=(abs(fft(w.*y))); max(fa), w=tukeywin(length(x_)); [f,Lf]=fourier_axis(x_); plot(f,abs(fa)); [mfa,mdx] = max(fa); [a,b,c]=gaussfit3(f(mdx+(-1:1)),fa(mdx+(-1:1))); hold on; plot(f,a*normpdf(f,b,c)); a_ = pi*a*L/length(x_)

c=lsqnonlin(@(c) c(1)*lpdf(f(mdx+(-1:1)),c(2),c(3)) - fa(mdx+(-1:1)),[fa(mdx),f(mdx),0.01]); plot(f,c(1)*lpdf(f,c(2),c(3)))


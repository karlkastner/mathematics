 y=[];
n=1e6;
 p=0.01;
L=2;
 x= linspace(-L/2,L/2,n)';
 y(:,1)=rectwin(x,0,p*0.75);
 y(:,2)=trapwin(x,0,p);
L_ = 1/sqrt(2)*0.75/2;
 s=p*sqrt(-0.5/log(0.5))*L_;
 y(:,3)=normpdf(x,0,s);
 y_=tukeywin(p*n/2+2);
 y(n/2-p*n/4:n/2+p*n/4+1,4)=y_;
S=abs(fft(y)).^2;
fx = fourier_axis(x);
clf
subplot(2,1,1)
 plot(fx,S./max(S));
xlim(3*[-100,100])
%plot(x,y./max(y)); xlim([-0.1,0.1])
if (0)
z = randn(n,1);
z = conv(z,y(:,1),'same');
S = periodogram(z,L);
hold on;
S = ifftshift(trifilt1(fftshift(S),101));
plot(fx,S./sum(S))
end
subplot(2,1,2)
plot(x,y./max(y))

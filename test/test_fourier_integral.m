% Wed 26 Jan 10:34:23 CET 2022

n = 1024;
L = 1;
%2*pi;
x = innerspace(-L/2,L/2,n)';
e = randn(n,1);
e = e-mean(e);
y = e;
z = e;
dx = (x(2)-x(1));

for idx=2:3
e(:,idx) = (cumsum(e(:,idx-1))-0.5*e(:,idx-1))*sqrt(dx);
%e(:,idx) = (cumsum(e(:,idx-1))-0*e(:,idx-1))*sqrt(dx);
end
% y = 1/6*abs(x).^3;
%y=y-mean(y);

fx=fourier_axis(L,n);
D = 2i*pi*fx;
I = 1./D;
I(1)=0;


if (0)
for idx=1:3
	subplot(2,3,idx);
	plot(x,y);
	if (idx<3)
	y = ifft(D.*fft(y));
	y = real(y);
	end
end
end
clf
for idx=1:3
	subplot(2,3,7-+idx);
	%plot(x,[y-mean(y),e(:,idx)-mean(e(:,idx)),z]);
	plot(x,[y-mean(y),e(:,idx)-mean(e(:,idx)),z-mean(z)]);
	E(:,idx) = y;
	%y = ifft(I.*fft(y));
	y = brownian_noise1d(L,n,y);
	y = real(y);
	z = brownian_noise_1d_acf(L,n,z);
end


%dy_dx = ifft(D.*fft(y));
%d2y = ifft(D.*fft(dy_dx));
%Y=ifft(I.*fft(y));
%y_ = ifft(I.*fft(dy_dx));
%plot(x,[y,y_,real(dy_dx),NaN*real(d2y),Y,1/3*x.^3]) 

% n=1e3; e = brownian_noise_1d_acf(1,[n,1e4]); x = linspace(0,1,n)'; plot([std(e,[],2),sqrt(x)])


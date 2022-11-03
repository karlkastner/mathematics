% Tue  3 May 20:19:15 CEST 2022
if (0)
n=1e4+1; x = zeros(n,1); x((end+1)/2) = 1; for idx=1:4; x = meanfilt1(x,round(n*1/5)); end; y=[x,fourier_derivative2(x,1),fourier_derivative2(x,2), fourier_derivative2(x,3)]; clf; plot(y./rms(y),'.-'); y = [x cdiff(x) cdiff(cdiff(x))]; %hold on;plot(y./rms(y),'--')

end

% TODO, this should be fourier diffusion, not weight

L=1;
 n = 1e2;
 x = randn(n,1).^2;
x = 0.*x;
x(end/2) = 1;
% w = caesaro_weight(L,n);
x_ = linspace(-5,5,n)';
dx = x_(2)-x_(1);
s = dx*[0.5,1,2,3,4];
m = length(s);
for idx=1:length(s)
 %x(:,idx+1) = meanfilt1(x(:,idx),3);
% x(:,idx+1) = trifilt1(x(:,idx),3);
 p = normpdf(x_,0,s(idx));
 x(:,idx+1) = conv(x(:,1),p,'same');
end

y=fourier_derivative(x,L,2);
% z=ifft(w.^10.*fft(x));
% y=fourier_derivative(z,2);
% min(x);
% x=(real([y,z]));
% plot(x./rms(x))
subplot(2,1,1)
plot(x)
for idx=1:m+1
subplot(2,m+1,m+1+idx)
y = real(y);
plot(y(:,idx)./rms(y(:,idx)))
end

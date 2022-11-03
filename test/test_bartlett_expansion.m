% 2022-04-01 17:59:13.503138871 +0200

if (0)
n = 1e5;
y = randn(n,1);
L = 100;
m = 2;
y = bandpass1d_implicit(y,0.9,1);

S = periodogram_bartlett(y,L,m); %,ni,pwin,subtract_mean)
n_ = n/m;
%f = [f(1:n_/2); zeros(n/2-n_/2,1); f(end-n_/2+1:end);zeros(n/2-n_/2,1)];

if (1)
f = fft(S);
f = padd_(f,n-n_);
subplot(2,2,3)
plot(f)
S_ = ifft(f)*m;
else
f = fft(S,n);
S_ = ifft(f(1:2*n_),n)*m;
end
S_(:,2) = periodogram_bartlett(y,L,m,n);
%S_(:,3) = periodogram_bartlett(y,L,1,n);
%S_(:,3) = abs(fft(y)).^2;

subplot(2,2,1);
plot(S)
subplot(2,2,2);
plot(S_)
end

% (F padd y)^2 + (F padd z)^2 =? F padd F^-1 (|F y|^2 + |F z|^2)
n=1e3;
m = 1e1;
 x = randn(n,m);
flag = 0;
for idx=1:m
 x(:,idx) = bandpass1d_implicit(x(:,idx),0.9,1);
end
 f = fft(x,m*n);
 S = sum(abs(f).^2,2);

 % joint interpolation
 f = fft(x);
 S_ = sum(abs(f).^2,2);
 f_     = fft(S_);
 f_     = padd_(f_,(m-1)*n);
 S(:,2) = ifft(f_)*m;

% individual interpolation yields same result as joint interpolation
% padding has therfore to the pattern before computing the periodogram
f  = fft(x);
S_ = abs(f).^2;
f_     = fft(S_);
f_     = padd_(f_,(m-1)*n);
S(:,3) = sum(ifft(f_)*m,2);



% a "trick" is to interpolate the square root
% joint interpolation
 f = fft(x);
 S_ = sum(abs(f).^2,2);
 f_     = fft(sqrt(S_));
 f_     = padd_(f_,(m-1)*n);
 S(:,4) = ifft(f_).^2*m^2;
clf
plot(real(S))



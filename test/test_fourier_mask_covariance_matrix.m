% cov((F*msk*x)*conj(F*msk*x), (F*msk*x)*conj(F*msk*x))

F msk x = F msk F.' F.'^-1 x

if (0)
n=1e2;
m = 1e4;
 mm=11;
 x = randn(n,m);
 x(n+1:end,:) = 0;
 msk=ones(n,1);
 msk(mm+1:end)=0;
 S = abs(fft(msk.*x)).^2;
 S_=S(1:n/2,:);
 C = cov(S');
 subplot(2,2,1);
 imagesc(C);
 colorbar;
% subplot(2,2,2);
% imagesc(abs(fft2(y'*y)))
subplot(2,2,2);
w = abs(fft(msk)).^2;
plot([w(:)/w(1,1),C(:,1)/C(1,1)])
CC = zeros(n,n);

subplot(2,3,5)
fC = abs(fft2(C));
imagesc(fC)

subplot(2,3,6)
plot(diag(fC)/fC(1,1))

end
if (1)
% 2D
n = 4e1;
m = 1e3;
mm = 1e1;
mm = n/10;

x = randn(n,n,m);
msk = zeros(n,n);
msk(1:mm,1:mm) = 1;

S = abs(fft2(msk.*x)).^2;
S = reshape(S,n*n,m);
C = corr(S');
subplot(2,2,3);
imagesc(C)
subplot(2,2,4);
w = abs(fft2(msk)).^2;
plot([w(:)/w(1,1),C(:,1)/C(1,1)])

subplot(2,3,5)
fC = abs(fft2(C));
imagesc(log10(fC))

subplot(2,3,6)
plot(diag(fC)/fC(1,1))

end

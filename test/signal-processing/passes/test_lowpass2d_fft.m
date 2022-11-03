% 2021-07-12 00:13:18.778834439 +0200

if (1)
 rho=[0.5,0.25];
 rho = 0.75;
 n= 4*[100,100];
 a = deg2rad(-45);
% n= [10,4];
 x = randn(n);
% y = lp2d_fft(x,rho,a);
% y = bandpass2d_fft(y,rho,a);
order = [2,4,6];
clf
for idx=1:length(order)
 y = x;
subplot(2,3,idx)
if (1)
 y = bandpass2d_implicit(y,rho,a,order(idx));
else
 rho = 0.5;
 y = lowpass2d_implicit(y,rho,a,order(idx)); %,a,order(idx));
end
 f = fftshift(abs(fft2(y)));
 % f = lowpass2d_implicit(f,0.5);
 imagesc(f);
 colorbar();
 axis equal;
 axis tight;
 f_C{idx} = f/max(f(:));
end
subplot(2,2,3)
f = triu(f_C{1}) + tril(f_C{end});
%imagesc(log10(f_C{1}) - log10(f_C{2}))
imagesc(f);
colorbar
axis equal
else
 rho = 0.5;
 n = 1e3
 x = randn(n,1);
 y = lp1d_fft(x,rho);
 plot(abs(fft(y))) 
end

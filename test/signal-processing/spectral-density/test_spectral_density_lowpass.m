% 2021-11-30 17:47:07.397464184 +0100
 p=1;
 n=20;
 x=0:n-1;
 fx=fourier_axis(x);
 r=0.8;
 x = zeros(n,1);
 x(1) =1;
 for idx=1:p;
 x = lowpass1d_implicit(x,[0,r],1,true);
 end;
  S=real(fft(x));
 S(:,2) = spectral_density_lowpass_one_sided(0:n-1,r,p,n);
 S=S./S(1,:);

% two-sided
x2 = zeros(n,1);
x2(1) =1;
for idx=1:p;
	x2 = lowpass1d_implicit(x2,[r,r],1,true);
end
S2=real(fft(x2));
S2(:,2) = spectral_density_lowpass(fx,r,p/2,1,'rho');
S2=S2./S2(1,:);


figure(1)
clf
subplot(2,3,1)
plot(x);

subplot(2,3,2)
plot(S)
ylim([0,1])

subplot(2,3,3)
plot(real(fft(S)));
hold on
plot(imag(fft(S)),'--');
%plot(x);

subplot(2,3,4)
plot([x2/x2(1)])
hold on
a=[r.^(0:19)' + r.^(20:-1:1)'];
plot(a/a(1))
plot(x/x(1))

subplot(2,3,5)
plot(S2)
ylim([0,1])

subplot(2,3,6);
plot(real(fft(S2)));
hold on
plot(imag(fft(S2)));

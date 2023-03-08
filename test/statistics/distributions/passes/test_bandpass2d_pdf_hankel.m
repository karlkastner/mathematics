%n=2.^(0:6);
%S=[];
L =50;
n = L^2;
fx = fourier_axis(n/L,n);
%x=linspace(0,4)';
a     = 0.5;
order = 1;
S = [];
t = [];
if (1)
tic()
S      = lowpass2d_pdf(fx,a,order);
t(1) = toc; tic;
[S(:,2),Rbar] = lowpass2d_pdf_hankel(fx,a,order);
t(2) = toc; tic;
[S2d,R2d]    = lowpass2d_discrete_pdf(L,n,1/a,order);
t(3) = toc; tic;
Rbar(:,2) = mean(R2d,2);
S(:,3) = S2d(:,1);
S(:,4) = real(ifft(mean(R2d,2)));
S=S./(S(1,:));
else
S = bandpass2d_pdf(x,a,order);
S(:,2) = bandpass2d_pdf_hankel(fx,a,order);
end

clf();
subplot(2,1,2)
plot(fx,S)
xlim([0,1/a])

subplot(2,2,2)
R = fft(S,[],1);
R = real(R);
plot(R./R(1,:));
subplot(2,2,3)
plot(Rbar./Rbar(1,:));
%xlim([0,1/a])
%for idx=1:length(n);
%	f=
%	S(:,idx) = f*w;
%end;
%S;
% resn = rms(S-S(:,end)), loglog(n,resn,'.')

rms(S-S(:,1))./rms(S)
rms(R-R(:,1))./rms(R)
t

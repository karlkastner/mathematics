% 2023-03-04 14:35:41.311550593 +0100
%n=2.^(0:6);
%S=[];
a     = 2;
m = 20;
L = m*2;
n = m^2;
fx = fourier_axis(L,n);
%x=linspace(0,4)';
order = 3;
S = [];
t = [];

tic()
S      = bandpass2d_pdf_exact(fx,a,order);
t(1) = toc; tic;
[S(:,2)] = bandpass2d_pdf_hankel(L,n,a,order);
t(2) = toc; tic;
[S2d,R2d]    = bandpass2d_discrete_pdf([L,L],[n,n],1/a,order);
t(3) = toc; tic;
%Rbar(:,2) = mean(R2d,2);
S(:,3) = S2d(:,1);
%S(:,4) = real(ifft(mean(R2d,2)));
%S=S./(S(1,:));
%S = bandpass2d_pdf(x,a,order);
%S(:,2) = bandpass2d_pdf_hankel(fx,a,order);

S = S./sum(S);

clf();
subplot(2,1,1)
plot(fx,S)
xlim([0,n/(2*L)])
rmse_S = rms(S(:,1:2)-S(:,1))./rms(S(:,1))
vline(a)

if (0)
subplot(2,2,2)
R = fft(S,[],1);
R = real(R);
plot(R./R(1,:));
subplot(2,2,3)
plot(Rbar./Rbar(1,:));
rms(R-R(:,1))./rms(R)
end
%xlim([0,1/a])
%for idx=1:length(n);
%	f=
%	S(:,idx) = f*w;
%end;
%S;
% resn = rms(S-S(:,end)), loglog(n,resn,'.')


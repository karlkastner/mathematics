% 2023-03-04 14:35:41.311550593 +0100
%n=2.^(0:6);
%S=[];

a   = 2;
m = 20;
L = m/2;
n = m^2; %L^2;
fx = fourier_axis(L,n);
%x=linspace(0,4)';
order = 1;
S = [];
t = [];

tic()
S      = lowpass2d_pdf(fx,a,order);
t(1) = toc; tic;
[S(:,2),Rbar] = lowpass2d_pdf_hankel(L,n,a,order);
t(2) = toc; tic;
[S2d,R2d,x,y] = lowpass2d_discrete_pdf([L,L],[n,n],1/a,order);

t(3) = toc; tic;
Rbar(:,2) = mean(R2d,2);
S(:,3) = S2d(:,1);

R = [];

%S(:,4) = real(ifft(mean(R2d,2)));
%S=S./(S(1,:));
r = hypot(cvec(x),rvec(y));
R(:,3) = R2d(:,1);
for idx=1:2
	fdx = (x>=0);
max(x)
max(r(:))
	S2d = interp1(x(fdx),S(fdx,idx),r,'spline',0);
	R2d = ifft2(S2d);
	R(:,idx) = R2d(:,1);
end

clf();
subplot(2,2,1)
plot(fx,S(:,1))
hold on
plot(fx,S(:,2:3),'--')
legend('exact','hankel','discrete');
xlim([0,1/a])
ylabel('S')
grid on
err_S = rms(S(:,1:2)-S(:,3))

subplot(2,2,2)
%R = fft(S,[],1);
R = real(R);
R = R./R(1,:);
plot(x,R(:,1))
hold on
plot(x,R(:,2:3),'--')
grid on
xlim([0,0.5*a]);
if (0)
subplot(2,2,3)
plot(Rbar./Rbar(1,:));
%xlim([0,1/a])
%for idx=1:length(n);
%	f=
%	S(:,idx) = f*w;
%end;
%S;
% resn = rms(S-S(:,end)), loglog(n,resn,'.')
end

rms(S-S(:,1))./rms(S)
rms(R-R(:,1))./rms(R)
t

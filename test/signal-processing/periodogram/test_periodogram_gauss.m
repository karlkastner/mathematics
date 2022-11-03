% Tue 17 May 15:13:49 CEST 2022
if (0)
L = 100; % 8e3
m = 100; %537;
f = 1; %0.01;
n = 50*L;
else
	L = 4e3;
	n = 4*L;
	f = 0.01;
	m = 10;	
end
x = linspace(-L/2,L/2,n)';
fx = fourier_axis(x);
y = sin(2*pi*f*x);
S      = periodogram(y,L);
S(:,2) = periodogram_gauss(y,L,m);
S(:,3) = periodogram_rectwin(y,L,m);
%S(:,4) = periodogram_bartlett(y,L,m,n);

Lw = m/L;
dx=x(2)-x(1);

fcut = rectwin_cutoff_frequency(Lw);
Lw_gauss = fcut2Lw_gausswin(fcut);
if (0)
	gw = normpdf(x,0,Lw_gauss);
else
	gw = normpdf(x,0,fcut2Lw_gausswin(rectwin_cutoff_frequency(1))*m*dx); % Lw_gauss
	rw = rectwin(x,0,m*dx);
%	rw = rectwin(fftshift(fx),0,m/L);
	%rw(:,2) = rectwin(x,0,Lw);
%clf
%plot(rw)
%pause
end
%id = (-n/2:n/2-1);
%rw = (abs(id)<=0.5*m);

S_ = [];
S_(:,2) = ifftshift(conv(fftshift(S(:,1)),gw,'same'));
S_(:,3) = ifftshift(conv(fftshift(S(:,1)),rw,'same'));
%S_(:,3) = ifftshift(meanfilt1(fftshift(S(:,1)),m));

%fftshift(meanfilt1(fftshift(S(:,1)),m));
S=S./sum(S);
mS = max(S(:,3));
S = S/mS;
S_=S_./sum(S_);
S_ = S_/mS;
%S=S./max(S);
fx=fourier_axis(L,n);
clf
subplot(2,2,1)
plot(fx,S);
hold on
%set(gca,'colororderindex',1)
plot(fx,S_,'--');
xlim([0.5*f,f*1.5])
legend('p','gws','rws','','gw','rw')
%rms(S(:,3)-S(:,4))/rms(S(:,3))
subplot(2,2,2)
plot(S);
hold on
plot(S_,'--')
legend('p','gws','rws','','gw','rw')

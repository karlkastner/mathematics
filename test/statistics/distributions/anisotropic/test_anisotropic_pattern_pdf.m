% 2023-03-13 22:15:51.218584715 +0100
k = 16;

dx = 1/k;

f0 = 1;
L = k*[1,1];
n = L(1)/dx*[1,1];
df = 1./L;

%sxy =0.25*[1,1];
sxy =0.125*[1,1];

m    = 25;
if (~exist('Shat','var'))
Shat = 0;
rng(0)
for idx=1:m
	idx
	[b,xy,S,f] = anisotropic_pattern(L,n,f0,sxy,'exact');
	Shat_ = abs(fft2(b-mean(b,'all'))).^2;
	Shat = Shat + Shat_;
end
Shat = 2*Shat/(sum(Shat,'all')*df(1)*df(2));
Rhat = real(ifft2(Shat));
Rhat = Rhat/Rhat(1);
end

Shatx=sum(Shat,1)'*df(2);
Shaty=sum(Shat,2)*df(1);
Rhatx=mean(Rhat,1)';
Rhaty=mean(Rhat,2);

S  = anisotropic_pattern_pdf(L,n,f0,sxy);
R  = anisotropic_pattern_acf(L,n,f0,sxy);

Sx=sum(S,1)'*df(2);
Sy=sum(S,2)*df(1);
Rx=mean(R,1)'; %'*df(2);
Ry=mean(R,2);%*df(1);
Rx = Rx/Rx(1);
Ry = Ry/Ry(1);
Rhatx = Rhatx/Rhatx(1);
Rhaty = Rhaty/Rhaty(1);

Sx_ = phase_drift_pdf(f.x,f0,sxy(1));
Rx_ = phase_drift_acf(f.x,f0,sxy(1));

% where does the scale factor come from?
scale__ = 4.5;
Sy_ = phase_drift_parallel_pdf(f.y,scale__*sxy(2));
Ry_ = phase_drift_parallel_acf(f.y,scale__*sxy(2));

if (0)
Sx__ = S(1,:)'*2/(sum(S(1,:))*df(1));
Shatx__ = Shat(1,:)'*2/(sum(Shat(1,:))*df(1));
Sy__ = S(:,1)*2/(sum(S(:,1))*df(1));
Shaty__ = Shat(:,1)*2/(sum(Shat(:,1))*df(1));
else
Sx__ = [];
Sy__ = [];
Shatx__ = [];
Shaty__=[]
end

subplot(2,4,1)
%imagesc(S)
%subplot(2,4,2)
plot(f.x,[Sx,Sx_,Shatx])
xlim([0,0.5*n(1)/L(1)]);
legend('2d analytic','1d analytic','2d empir')
title('$\bar S_x$','interpreter','latex')
xlim([0,f0*2.5])

subplot(2,4,2);
plot(f.y,[Sy,Sy_,Shaty])
xlim([0,0.5*n(1)/L(1)]);
legend('2d analytic','1d analytic','2d empir')
xlim([0,f0*2.5])
title('$\bar S_y$','interpreter','latex')

if (0)
subplot(2,4,3);
plot(f.x,[Sx__,Sx_,Shatx__])
title('$S(x,0)$','interpreter','latex')
xlim([0,f0*2.5])

subplot(2,4,4);
plot(f.x,[Sy__,Sy_,Shaty__])
title('$S(y,0)$','interpreter','latex')
xlim([0,f0*2.5])
end

%subplot(2,4,5)
%imagesc(Shat)

subplot(2,4,5)
plot(xy.x,[Rx,Rx_,Rhatx])

subplot(2,4,6);
plot(xy.y,[Ry,Ry_,Rhaty])

subplot(2,4,7)
plot(xy.x,[R(1,:)',Rx_,Rhat(1,:)'])

subplot(2,4,8)
plot(xy.x,[R(:,1),Ry_,Rhat(:,1)])


figure(1)
clf
p = 2;
m = 200;
normalize = true;
f0=2;
 n=1e5;
 L = 100*[1,10,100];

% note that fmax for bartlett is 1/2*(n/m)/L !

%rho = 0.9;
S = [];
S_ = [];
Sb = [];
Sf = [];
fx = [];
rho = filter_f0_to_rho(f0,L,n)
%rho = rho(1).^(L)
%pause
if (1)
 for idx=1:length(L);
%r      = rho./(1-2*rho+rho.*rho);

 x=linspace(0,L(idx),n)';
 fx(:,idx)=fourier_axis(x);
 S(:,idx)=spectral_density_bp(fx(:,idx),f0,L(idx),length(x),p,'f',normalize);
% S_(:,idx)=spectral_density_bp(fx(:,idx),rho(idx),L(idx),length(x),p,'rho',normalize);
 S_(:,idx)=spectral_density_bp_approx(fx(:,idx),f0,L(idx),length(x),p,'f',normalize);
 y = randn(n,1);
 y = bandpass1d_implicit(y,rho(idx),p,true);
 Sf(:,idx) = meanfilt1(periodogram(y,L(idx)),m);
 Sb(:,idx) = periodogram_bartlett(y,L(idx),m,n);

% sum(S)/L(idx),
%max(fx)
 end;
subplot(2,3,1)
plot(fx,S);
xlim([0,4*f0]);
%hold on
%hold on
%plot(fx,S_,'--','linewidth',2)
subplot(2,3,2)
plot(fx,Sf);
 xlim([0,4*f0]);
subplot(2,3,3)
%plot(fx,Sb);
plot(fx,Sb);
 xlim([0,4*f0]);


end

if (0)
f0=2;
 N=[1e2,1e3,1e4];
 m =20;
 L = 10;
 L = 1;
 f0 = 10;
rho = filter_f0_to_rho(f0,L,N)
 for idx=1:length(N);
 n=N(idx);
 x=linspace(0,L,n)';
 fx=fourier_axis(x);
 S=spectral_density_bp(fx,f0,L,length(x),p,'f',normalize);
 y = randn(n,1);
 y = bandpass1d_implicit(y,rho(idx));
 Sb = periodogram_bartlett(y,L,m,n);
 sum(S)/L
 max(fx)
subplot(2,2,3)
 plot(fx,S);
 hold on;
 xlim([0,4*f0]);
subplot(2,2,4)
 plot(fx,Sb);
 hold on;
 end;
 xlim([0,4*f0]);
end

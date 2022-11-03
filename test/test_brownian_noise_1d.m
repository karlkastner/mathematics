% Tue 25 Jan 22:08:24 CET 2022
% The integral has to be periodic, so noise can only be integrated to a brownian
% bridge

L = 1;
 n = [1e2,1e4];

 e0=randn(n);
e0_ = e0;
%e0 =e0-mean(e0);
e0(end,:) = -sum(e0(1:end-1,:));

 % integrate in real space
dx=L/n(1);
 e_ = sqrt(dx)*(cumsum(e0)-0.5*e0);
 e__ = sqrt(dx)*(cumsum(e0_)-0.5*e0_);
 
 [e,Sb] = brownian_noise1d(L,n,e0);

%e = e+mean(e_);
%e = e+e_(:,end);

 clf;
subplot(2,2,1)
 plot([mean(e_,2),mean(e,2)]);
 hold on;

subplot(2,2,2)
 plot([std(e__,[],2),std(e_,[],2),std(e,[],2)])

subplot(2,2,3);
 fx=fourier_axis(L,n(1));
 S = mean(periodogram(e_,L),2);
 S(:,2) = mean(periodogram(e,L),2);
% S(:,2) = 1./(pi*fx.^2);
% S(:,3) = abs(Sb).^2;
%S(~isfinite(S)) = 0;
%S = S./sum(S);
 plot(fx,S);
%[mean(S,2),1./(pi*fx.^2)]);

sum(S) 

subplot(2,2,4)
plot([e_(:,end),e(:,end)])

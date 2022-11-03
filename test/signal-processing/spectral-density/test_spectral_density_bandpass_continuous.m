% 2021-11-23 17:54:36.646762263 +0100

NN  = [10,30,100,300,1e3]+1;
LL = round([10,30,100,300,1e3]/2);
%LL =  sqrt(N);
f0 = 1;
p  = 1;
figure(1)
clf
rmse = [];
for idx=1:length(N)
for jdx=1:length(LL)
	n = NN(idx);
	L = LL(jdx);
	dx = L/n;
	%dx = 1/sqrt(n);
%dx*n;
	df = 1/L;
	%fmax = 1/2*1/dx;
	%x = L*(0:n-1)/n;
	fx = fourier_axis(L,n);
	fdx = fx>0;
	%f = (0:df:fmax)';
	S = spectral_density_bandpass_discrete(fx,f0,p,dx,[],'f')
	S(:,2) = spectral_density_bandpass_continuous(fx,f0,p);
	S = S./(sum(S(fdx,:))*df);
	subplot(nN,nN,(idx-1)*nN+jdx);
	cla;
	plot(fx(fdx),S(fdx,:));
	rmse(idx,jdx) = rms(diff(S(fdx,:),[],2));
	title(sprintf('n %d L %g %g%%',[n,L,1e2*rmse(idx,jdx)]));
end
end
figure(2);
clf
subplot(2,2,1)
loglog(NN,(rmse),'o-')
legend(num2str(cvec(LL)))
subplot(2,2,2)
loglog(LL,(rmse)','o-')
legend(num2str(cvec(NN)))

subplot(2,2,3)
dx  =rvec(LL)./cvec(NN);
loglog(dx,rmse,'o-')
dx_ = logspace(-3,2);
hold on
plot(dx_,dx_.^2)
title('rmse, quadratic dx')

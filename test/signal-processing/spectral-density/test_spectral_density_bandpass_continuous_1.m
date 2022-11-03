% Tue 23 Aug 13:57:17 CEST 2022

fx = linspace(0,100,1e5)';
p  = [1,2];
fc = [1,2];
S = [];
Sn = [];
k = 0;
for idx=1:length(p)
for jdx=1:length(fc)
	k = k+1;
	S(:,k)  = spectral_density_bandpass_continuous(fx,fc(jdx),p(idx),false);
	Sn(:,k) = spectral_density_bandpass_continuous(fx,fc(jdx),p(idx),true);
end
end

IS = cumsum(mid(Sn).*diff(fx));

figure(1);
clf
subplot(2,2,1)
%loglog(fx,S);
plot(fx,S);
xlim([0,10])
subplot(2,2,2)
plot(mid(fx),IS)



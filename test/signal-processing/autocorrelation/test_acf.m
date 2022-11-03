% 2021-09-03 17:34:35.444456536 +0200

L = 100;
Lp = 2;
m=8;
x = linspace(0,L,1e5)';
clf;
w = [1,rand()];
for idx=1:m
	k = 2*pi*[1,2];
	phi = [0,2*pi*rand(1,1)]; %(idx-1)/4*pi];
	subplot(2,m,idx)
	y = sum(w.*cos(x*k+phi),2);
	% + sin(2*x+[(0:4)*pi/4]);
	plot(x,y,'-')
	ylim([-2,2])
	xlim([0,Lp])
	subplot(2,m,m+idx);
	a = autocorr_fft(y,[],true);
	%a = -sum(1./k.*sin(x*k+phi),2);
	a(:,2)= sum(w.^2.*cos(k.*x),2)/sum(w.^2);
	plot(x,a);
	%/sum(w.^2);
	%hold on;
	%plot(x,a);
	xlim([0,Lp])
	ylim([-1,1])
end


% Wed 22 Sep 13:26:33 CEST 2021

n=1e3;
nf =11;
fmin = 1/200;
fmax = 11/200;
L=2e3;
x=linspace(0,L,n)';
fx=fourier_axis(x);

nn = 1e3;
P = [0.5,0.2,0.1,0.05,0.02,0.01];
if (1)
disp('false positive')
fdx = find(fx>1/200 & fx<11/200);
pp = zeros(4,length(P));
for idx=1:nn
	y = 0*cos(2*pi*x/100) + 1*randn(n,1);
	y = y/std(y);
	y= y-mean(y);
	S=periodogram(y,L);
	Se = spectral_density_iid(L,n);
	Seb = periodogram_bartlett(y,L,nf,n);
%	Seb = nf/(nf-1)*(Seb - S/nf);
	Sme = median_periodogram(S,nf);
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'exact',Se);
	pp(1,:) = pp(1,:) + (p<P)/nn;
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'inclusive',nf);
	pp(2,:) = pp(2,:) + (p<P)/nn;
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'exclusive',nf);
	pp(3,:) = pp(3,:) + (p<P)/nn;
	p = periodogram_p_value(max(S(fdx)),Se(1),1);
	pp(4,:) = pp(4,:) + (p<P)/nn;
end
[P; pp]
clf()
plot(P,pp)
%[P',ppa]

end
if (1)
disp('false negative')
pp = zeros(4,length(P));
for idx=1:nn
	y   = 1*cos(2*pi*x/100) + 0.0*randn(n,1);
	y   = y/std(y);
	y   = y-mean(y);
	S   = periodogram(y,L);
	Se = spectral_density_iid(L,n);
	Seb = periodogram_bartlett(y,L,nf,n);
	% Sme = median_periodogram(S,nf);
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'exact',Se);
	pp(1,:) = pp(1,:) + (p>P)/nn;
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'inclusive',nf);
	pp(2,:) = pp(2,:) + (p>P)/nn;
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,S,fmin,fmax,'exclusive',nf);
	pp(3,:) = pp(3,:) + (p>P)/nn;
	p = periodogram_p_value(max(S(fdx)),Se(1),1);
	pp(4,:) = pp(4,:) + (p>P)/nn;
end
[P;pp]

end

if (1)
disp('power');
% TODO, check that the correct peak is identified (!)
% note : exclusive detects peak more often, but also other peaks
%se = [0.01,0.1,1,10,100];
se = [1,10:10:90]/10;

P = 0.95;
pp = zeros(4,length(se));
for jdx=1:length(se)
for idx=1:nn
	y   = 1*cos(2*pi*x/100) + se(jdx)*randn(n,1);
	y   = y/std(y);
	y   = y-mean(y);
	Shat   = periodogram(y,L);
	Se  = 2*L/n*ones(n,1);
	Seb = periodogram_bartlett(y,L,nf,n);
	% Sme = median_periodogram(S,nf);
	[p,ratio,maxShat,mdx,fdx] = periodogram_test(fx,Shat,fmin,fmax,'exact',Se);
	pp(1,jdx) = pp(1,jdx) + (p<P)/nn;
	pp_(1,jdx,idx) = p;
	[p,ratio,maxShat,mdx,fdx,S_in] = periodogram_test(fx,Shat,fmin,fmax,'inclusive',nf);
	pp(2,jdx) = pp(2,jdx) + (p<P)/nn;
	pp_(2,jdx,idx) = p;
	[p,ratio,maxShat,mdx,fdx,S_ex] = periodogram_test(fx,Shat,fmin,fmax,'exclusive',nf);
	pp(3,jdx) = pp(3,jdx) + (p<P)/nn;
	pp_(3,jdx,idx) = p;
	p = periodogram_p_value(max(Shat(fdx)),Se(1),1);
	pp(4,jdx) = pp(4,jdx) + (p<P)/nn;
	pp_(4,jdx,idx) = p;
end % for idx
figure(2)
subplot(2,5,jdx)
cla
fdx = fx>0;
n_ = sum(fdx);
a = 1;                                                                          
b = (nf-1);
p = (1:n_)/(n_+1);                                               
qbeta = nf*betainv(p,a,b); 
plot(qbeta,sort(Shat(fdx)./S_in(fdx)));
hold on
plot([0,5],[0,5])

figure(3)
subplot(2,5,jdx)
%plot([Shat, S,S_])
plot(fx(fdx),Shat(fdx)./[Se(fdx), S_in(fdx),S_ex(fdx)])
title(se(jdx))
%xlim([0,40])
%pause
end % for jdx
[se
median(pp_,3)]
[se;
 pp]
end
% x = mean(randn(1e6,11,2).^2,3); mean(mean(x,2)), std(mean(x,2)), mean(median(x,2)),std(median(x,2)) 
%n=1e4; p=0.95; c = mean(randn(1e4,2).^2,2); mean(c>1/2*chi2inv(p,2)), x = randn(n,1); x=x/std(x); c=abs(fft(x)).^2/n; mean(c>1/2*chi2inv(p,2))  


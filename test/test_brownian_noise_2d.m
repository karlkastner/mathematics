% Thu 27 Jan 09:24:30 CET 2022

% keyword : brownian noise, surface sheet

%n=[1000,1e3]/25+1;
n = [1,1]*100;
m = n(1);
L = 1;

if (0)
	fx = fourier_axis(2,n(1)+m);
	x=linspace(0,1*(1+m/n(1)),n(1)+m)';
	dx=x(2)-x(1);
	e = randn(n);
	e = [e;
	     -e];
	     % -mid(e), better than e(1)
	%     -ones(m,1)*sum(e,1)/m];
	sS = 1./fx;
	sS(1)=0;
	B=sqrt(n(1))/(2*pi)*ifft(sS.*fft(e));
	B=B-B(1,:);
	subplot(2,2,1)
	plot([std(B*sqrt(n(1)/(n(1)-1)).^0,[],2),sqrt(x)])
end

figure(1)
clf
figure(2)
clf
fx = fourier_axis(L,n(1)+m);
e = randn(n(1));
e = [e, -e];
e = [e; -e];
%e = randn(2*n(1));
for idx=1:3

switch (idx)
case {1}
	% 2D
	df = 1/L;
	fr = hypot(fx,fx');
%	fr = abs(fx) + abs(fx');
%	fr_ = (abs(fx).^2 + abs(fx').^2).^(1/2);
%	fr = (fr).^1.9;
%        fr = fr+0.5*df;
%	fr = fr.^1.4;
%	fr = fx .* (fx');
%	fr = abs(fr);
case {2}
	% along x
	fr = hypot(fx,0.*fx');
case {3}
	% along y
	fr = hypot(0.*fx,fx');
end

sS = 1./fr;
sS(~isfinite(sS))=0;

figure(1)
subplot(2,3,idx)
imagesc(log10(fftshift(sS)))

subplot(2,3,3+idx)
switch (idx)
case {1}
	R = ifft2(sS);
	plot(R(1,:)/R(1,1))
	hold on
	plot(R(:,1)/R(1,1))
%	imagesc(real(R))
case {2}
	R = ifft(sS(:,1).^2);
	R(:,2) = ifft(sS(1,:).^2);
	plot(R/R(1));
case {3}
	R = ifft(sS(:,1).^2);
	R(:,2) = ifft(sS(1,:).^2);
	plot(R/R(1));
end

figure(2)
switch (idx)
case {1}
if (0)
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n*2,[L,L],2,true);
	%[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n*2,[L,L],2,true);
	D2 = D2x + D2y;
tic
	B = (-D2)^(-0.5) * flat(e);
toc
%	B = (-D2) \ flat(e)
%B
%pause
	B = reshape(B,size(e));
else
	B = ifft2(sS.*fft2(e));
end
%	id = (0:2*n-1);
%	r = hypot(id,id');
%	B = B - B(1,1); %end/2,end/2);
%	[r,sdx] = sort(flat(r));
%	B_ = B(sdx);
	%B = B - B(end/2,:);
	%B = B - B(:,end/2);
case {2}
	B = ifft(sS.*fft(e));
	B = B - B(end/2,:);
case {3}
	B = ifft(sS.*fft(e,[],2),[],2);
	B = B - B(:,end/2);
end
B = real(B);

subplot(2,3,idx)
if (1)%2 == idx || 3 == idx)
%sS(1,1) = 0;
	%plot([std(B,[],1)', std(B,[],2)])
	%plot([var(B,[],1)', var(B,[],2)])
if (1)
	a=autocorr_fft(B);
	plot(mean(a,2));
	hold on
	a=autocorr_fft(B');
	plot(mean(a,2));
end

else
if(0)
	a=autocorr_fft(B);
	plot(mean(a,2));
	hold on
	a=autocorr_fft(B');
	plot(mean(a,2));
end
%ACCUMARRAY(IND,DATA,SZ,FUN)
if (1)
%-> this is clearly wrong, the variance does not increase linearly!

%	B_ = B - B(:,1);
%	y = mean(B_.^2);
%	plot(y)
%	hold on
%	B_ = B - B(1,:);
%	y = mean(B_.^2,2);
%	plot(y)
end
if (0)
cla
%imagesc(log10(fftshift(sS)))
plot(sS(1,:))
hold on
plot(sS(:,1))
plot(sS(:,end/2))
end
	
%	B_ = accumarray(round(flat(r))+1,flat(B),[],@(x) mean(x.^2));
%	plot(B_);
%meanfilt1(B_,200));
%	plot(B(1,:))
%	hold on
%	plot(B(:,1))
%	plot(r(:),B(:),'.')
end
subplot(2,3,idx+3)
%imagesc(B)
plot(sort(flat(B)))

figure(3)
subplot(2,3,idx)
cla
	s2 = [];
	B = B/std(B(:));
idx
	for jdx=1:size(B,2)-1
		s2(jdx,1) = mean(mean((B(:,jdx)-B(:,1)).^2));
		s2(jdx,2) = mean(mean((B(jdx,:)-B(1,:)).^2));
%		s2(jdx,1) = mean(mean((B(:,jdx+1:end)-B(:,1:end-jdx)).^2));
%		s2(jdx,2) = mean(mean((B(jdx+1:end,:)-B(1:end-jdx,:)).^2));
	end
%	clf
	plot(s2)
%	pause

end

if (0)
	L=10;
 fmax=10;
 syms i o t f;
 I=int(exp(1i*2*pi*f*t)*1/(2*pi*f)^2,f), f=(subs(subs(I,f,fmax))) + (subs(subs(I,f,-fmax)));
 f = matlabFunction(f);
 plot(imag(f(linspace(-L,L))))
end


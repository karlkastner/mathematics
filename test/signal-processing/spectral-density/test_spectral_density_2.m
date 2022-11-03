% Fri  8 Apr 11:49:15 CEST 2022

% results :
%	Sc/lambda = Sc*fc is larger in 2d than in 1d
%		theoretical limit for Sc = L1*L2, the peak in 2d is more pronounced
%		test comprises also of n1*n2 bins
%	Sc limit has more complicated expression, when domain is not square or
%		pattern not align with coordinate axes
%	via truncation of R, determination of m is difficult
%		
if (0)
img = imread('img/tiger.jpg');
x = [275 1168]
y = [103 287]

img = mean(img,3);

%round(limits(b))
%close all; imagesc(img)
%nc=length(c); S=periodogram_bartlett(c-mean(c),nc,7,nc); x=1:nc; fx=fourier_axis(x); plot(fx,S); [Sc,l] = max(S), xlim([0,3.5*fx(l)]);  Sc*fx(l) 

% rotate and crop pattern
alpha = -12;
crop = [376   271   797   339];
img_ = imcrop(imrotate(img,alpha),crop);
if (0)
imagesc(img_)
end
r0 = 10;
n = size(img_);
n(1)=n(2);
x = repmat(1:n(2),n(1),1);
phi = 0;
%phi = 2*pi*randn(n(1),1);
%phi = cumsum(phi)/20;
%rho = 0.1;
%phi = 2*pi*0.4*brownian_noise_2d_fft(n,n);
%phi =0.03*brownian_noise_2d_fft(n,n,phi);
%phi = filter(1-rho,[1,rho],phi);
%phi = filter(1-rho,[1,rho],phi')';
img_ = cos(2*pi*x/75 + phi);
f0 = 0.004;
[S,R,fr,fx,fy]=spectral_density_estimate_2d(img_,r0,f0);
clf;
imagesc(fftshift(S));
[Sc,mdx]=max(S(:))
fc = fr(mdx);
%Sc_=Sc/sqrt(prod(size(S)))
regularity = sqrt(Sc*fc^2)
end

C = {
'tiger.jpg',	12, [376   271   797   339],	0.004
'zebra.jpg',	-10, [1161.5 1016.5 1490 591],	0.004
'leopard-1.jpg', 0, [516.5 309.5 526 258],	0.01
'leopard-2.jpg', 7, [2331.5 1301.5 1194 680],	0.005
'fish-plectorhinchus-chaetodonoides.jpg', 0, [303.5 240.5 329 201], 0.01
'fish-holacanthus-ciliaris.jpg', -30, [963.5 1144.5 703 424], 0.02
}
for idx=1:size(C,1)
	img = imread(['img/',C{idx,1}]);
	img = imrotate(img,-C{idx,2});
	if(isempty(C{idx,3}))
		clf
		imcrop(img)
	else
		img = imcrop(img,C{idx,3});
	end
	figure(1)
	subplot(2,3,idx);
	imagesc(img)
	img = double(img);
	img = mean(img,3);
	img = img-mean(img(:));
%	Shat = abs(fft2(img)).^2;
%	Shat = fftshift(Shat);
%	nf = 5;
%	Shat = trifilt2(Shat,nf);
	r0 = 5;
	r0 = 50;
	r0 = min(size(img))/5
	r0 = 5;
	method = 'tri';
	r0 = 9;
	method = 'rect';
	%f0 = 0.004;
	f0 = C{idx,4};
	[S,R,fr,fx,fy]=spectral_density_estimate_2d(img,r0,f0,method);
	
	figure(2);
	subplot(2,3,idx);
	%imagesc(log10(Shat))
	imagesc(fftshift(fy),fftshift(fx),fftshift(S))
	axis equal
	axis(0.2*[-1,1,-1,1])
	[Sc(idx,1),mdx] = max(S(:));
	fc(idx,1) = fr(mdx);
end
Sc
fc
reg = sqrt(Sc).*fc;
figure(3);
bar(reg);
set(gca,'xtick',1:length(reg),'xticklabel',C(:,1))


% Mon 13 Mar 11:18:34 CET 2023

dx  = 0.1;
L   = 10*[1,1];
n   = L/dx;
Scx = 1;
Scy = Scx;
fc  = 1;

	[f0,sx] = phase_drift_pdf_mode2par(fc,Scx);
	sy      = phase_drift_parallel_pdf_max2par(Scy);
sx
sy
sy = sx;
m = 100;
if (~exist('Shat','var'))
Shat = 0;
for idx=1:m
	rng(idx)
	[be,~,~,f] = anisotropic_pattern(L,n,f0,[sx,sy],true);
	Shat = Shat + abs(fft2(be-mean(be(:)))).^2;
end
end

Sx = mean(Shat)';
fdx = (f.x>=0);
Sx = Sx./(sum(mid(Sx(fdx)).*diff(f.x(fdx))));

Sy = mean(Shat,2);
fdx = (f.y>=0);
Sy = Sy./(sum(mid(Sy(fdx)).*diff(f.y(fdx))));


figure(1);
clf
subplot(2,4,1)
imagesc(be>mean(be(:)))
subplot(2,4,2)
imagesc((fftshift(Shat)))
subplot(2,4,3)
plot(Sx);
subplot(2,4,4)
plot(Sy);

[b,xy,S,f] = anisotropic_pattern(L,n,f0,[sx,sy*4],false);

Sx = mean(S.S2d)';
fdx = (f.x>=0);
Sx = Sx./(sum(mid(Sx(fdx)).*diff(f.x(fdx))));

Sy = mean(S.S2d,2);
fdx = (f.y>=0);
Sy = Sy./(sum(mid(Sy(fdx)).*diff(f.y(fdx))));

subplot(2,4,5)
imagesc(b>mean(b(:)))

subplot(2,4,6)
imagesc((fftshift(S.S2d)))

subplot(2,4,3)
hold on
plot(Sx);
subplot(2,4,4)
hold on
plot(Sy);



% Tue  6 Dec 10:04:56 CET 2022
figure(1);
clf
d=5;
a = 0:d:(360-d)
slope = NaN(length(a),1);
sS = NaN(length(a),1);
for idx=1:length(a)
idx
%try
	subplot(ceil(sqrt(360/d)),ceil(sqrt(360/d)),idx)
	[fx,fy]= fourier_axis_2d([1,1],100*[1,1]);
	Sx = spectral_density_bandpass_continuous(fx,10,2);
 Sy=normpdf(fy,0,2);
 S=Sx*Sy';
S = S+0.0001*randn(size(S))
 S=ifftshift(imrotate(fftshift(S),a(idx),'crop'));
sS(idx,1) = sum(S(:));
 [S,~,~,slope(idx,1)] = separate_isotropic_from_anisotropic_density(S,[1,1])
imagesc(fftshift(S.hat));
title([a(idx),slope(idx,1)])
axis square
%catch
%end

end
figure(2)
plot(a',[tand(a'),slope]);
hline(1);


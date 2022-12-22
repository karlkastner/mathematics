L = [1,1];
n = [10,10]*10;
[fx,fy,fr] = fourier_axis_2d(L,n);
fc = 10;
S = spectral_density_bandpass_continuous(fr,fc,2);
Shat = S;
[Sr, fri, count, A] = periodogram_radial(Shat,L);

S_ = reshape(A*Sr.normalized,n);

subplot(2,2,1)
imagesc(S)
subplot(2,2,2)
imagesc(S_)
subplot(2,2,3)
plot(Sr.normalized);
subplot(2,2,4)
imagesc(S_/sum(S_(:))-S/sum(S(:)))



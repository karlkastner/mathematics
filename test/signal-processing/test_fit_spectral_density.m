% Mon  6 Feb 14:26:57 CET 2023
% Karl Kastner, Berlin

fc = 1;
L = 10;
n = L^2;
f = fourier_axis(L,n);
f = cvec(f);
p = 10;
w = 1;
%S = spectral_density_bandpass_continuous(f,fc,p);

s = 0.3;
S = phase_drift_pdf(f,fc,s);
Sfun = @phase_drift_pdf;

clf
plot(f,S)

m = 1e1;
% p  = 1000;
for idx=1:2
	w  = 1;
	% p can never be lower than 0.5 
	p0 = [0,0.5] + rand(1,2);
	p  = fit_spectral_density(f,S,w,Sfun,p0,'mise-cramer')
	nf = 3;
	Sp = Sfun(f,p(1),p(2));
	subplot(2,2,idx)
	plot([S,Sp])
	rms(S-Sp)
%
%	p  = fit_spectral_density(f,S,w,'bandpass-continuous',p0,'ls',nf)
%	S = S.*chi2rnd(m,n,1)/m;
end
% with noi


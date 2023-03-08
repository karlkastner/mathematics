% 2023-03-07 12:57:23.982954801 +0100
f0 = 2*pi;
a = 2.^(-2:2);


m = 20;
l0 = 6/f0;
dx = l0/m;
L = m*l0;
n = L/dx;
fx = fourier_axis(L,n)';
df = 1/L;

fc = []
Sc = []
for idx=1:length(a)
	[fc(idx,1),Sc(idx,1)] = bandpass2d_pdf_mode(f0,a(idx));
	S = bandpass2d_pdf(fx,f0,a(idx));
	S = 2*S/sum(S*df);
	subplot(2,3,idx)
	cla
	plot(fx,S);
	hold on
	plot(fc(idx,1),Sc(idx,1),'*')
end

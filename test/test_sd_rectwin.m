% Tue 17 May 15:03:27 CEST 2022

L = 20;
n = 1e3;
Lw = 0.2;
x = linspace(-L/2,L/2,n)';
fx = fourier_axis(x);
if (1)
	w = abs(x)<0.5*Lw;
	S = abs(fft(w-0*mean(w))).^2;
	S(:,2)   = sd_rectwin(L,n,Lw);
	fcut     = rectwin_cutoff_frequency(Lw);
	Lw_gauss = fcut2Lw_gausswin(fcut);
	%^S(:,3)=normpdf(fx,0,Lw_gauss);
	S(:,3) = sd_gausswin(L,n,Lw_gauss);
	S = S./max(S);

else
	w = normpdf(x,0,Lw);
	S = abs(fft(w-0*mean(w))).^2;
	S(:,2) = sd_gausswin(L,n,Lw);
	%S(:,2) = sd_gausswin(L,n,Lw);
	s = 1/(sqrt(8)*pi*Lw);       
	% log(0.5) = -0.5*x^2/(1/(sqrt(8)*pi*Lw)^2)
	% sqrt(-2*log(0.5))/ 
	S(:,3) = exp(-fx.^2/(2*s^2));
	S = S./max(S);
	%fcut = -(2*pi*Lw)*log(0.5);
end
plot(fx,S)
vline(fcut);
hline(0.5)
xlim([0,4*fcut])


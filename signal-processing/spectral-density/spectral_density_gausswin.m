function S = sd_gausswin(L,n,Lw)
	fx = fourier_axis(L,n);
	s = 1/(sqrt(8)*pi*Lw);
	S = normpdf(fx,0,s);
%	S = sinc(f.*Lw).^2
%	S(:,2) = (1 - x.^2/6).^2;
end


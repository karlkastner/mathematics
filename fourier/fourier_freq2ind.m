function f = fourier_freq2ind(f,n)
	fdx = f < 0;
	f(fdx) = f(fdx) + n+1;	
	f(~fdx) = f(~fdx) + 1;
end

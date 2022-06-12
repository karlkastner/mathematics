% Fri  8 Apr 11:40:40 CEST 2022
function [S,R,fr,fx1,fx2] = spectral_density_estimate_2d(img,r0,f0,method)
	Shat=abs(fft2(img-mean(img(:)))).^2;
	R = ifft2(Shat);
	n = size(R);
	L = n;
	fx1=fourier_axis(n(1),n(1));
	fx2=fourier_axis(n(2),n(2))';
	fr = hypot(fx1,fx2);
	x1 = fourier_axis(1,n(1));
	x2 = fourier_axis(1,n(2));
	r  = hypot(x1,x2');
	if (nargin()<4)
		method = 'trifilt';
	end

	switch(method)
	case {'acf'}
		R(r>r0) = 0;
		S = real(fft2(R));
	case {'rect'}
		%r0 = 5;
		S = Shat;
		S = fftshift(S);
		S = meanfilt2(S',round(r0))';
		S = ifftshift(S);
	case {'tri'}
		S = Shat;
		S = fftshift(S);
		S = trifilt2(S,round(2*r0));
		%S = trifilt1(S',round(2*r0))';
		S = ifftshift(S);
	end
	S(abs(fr)<f0) = 0; 
	S = (S);
	df = 1./L;
	S  = 2*S/sum(S(:)*df(1)*df(2));
end


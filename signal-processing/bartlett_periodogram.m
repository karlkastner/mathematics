% Fri  2 Jul 18:22:26 CEST 2021
% function S = bartlett_periodogram(fx,y,m,nn)
function S = bartlett_periodogram(fx,y,m,nn)
	if (isvector(y))
		y = cvec(y);
	end
	n = floor(size(y,1)/m);
	if (nargin()<4)
		nn = size(y,1);
	end
if (0)
	w = tukeywin(n);
	w = w/mean(w);
else
	w = 1;
end
	S = 0;
	for idx=1:m
		y_ = y((idx-1)*n+(1:n));
		y_ = y_-mean(y_);
		f  = fft(w.*y_,nn);
		S = S+ abs(f).^2/(n*mean(w.^2));	
	end
	S = S/m;

	df = fx(2)-fx(1);
	S = S/(2*df*sum(S));
end

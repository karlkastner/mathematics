% 2016-10-23 23:15:29.679880879 +0200
% Karl Kastner, Berlin
%% spectrogram
function [fy fx Tx Xi Yi] = spectrogram(X,Y)
	p = 0.5; % results in 1/4 of the range at each side
	n  = length(X);
	ni = 2.^ceil(log(n)/log(2));
	% interpolate to constant step width
	Xi = X(1) + (0:ni-1)/ni*(X(end)-X(1));
	Yi = cvec(interp1(X,Y,Xi));
	% apply spectral analysis piecewise
	% split the series into np pieces and estimate the periodogram individually
	w  = tukeywin(ni,p);
	w  = length(w)/sum(w)*w;
	fy = abs(fft(w.*Yi));
	[fx Tx mask] = fourier_axis(Xi);
	Tx  = Tx(mask);
	fy  = fy(mask);
	l = length(X); %X(end)-X(1);
	fy  = fy/l; %sqrt(l);
	fx  = fx(mask);
	
end


% Di 12. Jan 14:27:40 CET 2016
% Karl Kastner, Berlin
%% bartlet spectrogramm
%% TODO sliding window
function [Xi Yi Tx fy] = bartlett_spectrogram(X,Y,np)
	n  = length(X);
	ni = 2.^ceil(log(n)/log(2));
	% interpolate to constant step width
	Xi = X(1) + (0:ni-1)/ni*(X(end)-X(1));
	Yi = cvec(interp1(X,Y,Xi));
	w  = tukeywin(ni,1/8);
	Yi = w.*Yi;
	% apply spectral analysis piecewise
	% split the series into np pieces and estimate the periodogram individually
	fy = zeros(ni/np,1);
	for idx=1:np
		fy = fy + abs(fft(Yi(1+(idx-1)*ni/np:idx*ni/np))).^2;
	end
	fy = sqrt(fy);

	[fx Tx mask] = fourier_axis(Xi(1:ni/np));
	Tx  = Tx(mask);
	fy  = fy(mask);
end


% Di 12. Jan 14:52:17 CET 2016
% Karl Kastner, Berlin
%% welch spectrogram
function [fy fx fys Tx Xi Yi] = welch_spectrogram(X,Y,np)
	% order of detrending
	o = 0;

	p = 0.5; % results in 1/4 of the range at each side

	n  = length(X);
	% expand to powers of 2
	ni = 2.^ceil(log(n)/log(2));
	% interpolate to constant step width
	Xi = X(1) + (0:ni-1)/ni*(X(end)-X(1));
	Yi = cvec(interp1(X,Y,Xi));
	% apply spectral analysis piecewise
	% split the series into np pieces and estimate the periodogram individually
	n_ = round(ni/np);
	% the fourier transform of the shortened window is only of length n_,
	% but the expansion is made nonseless to ni, to spectrally interpolate
	% to the resolution of the unwindowed series
	ne = ni;
	fy = zeros(ne,1);
	% window for unperiodic data series
	w  = tukeywin(n_,p);
	w  = length(w)/sum(w)*w;
	for idx=1:2*np-1
		% let windows overlap by 50%, such that the tukey flaps are outside of the window range
		id = round(cvec(1+(idx-1)*ni/(2*np):(idx+1)*ni/(2*np)));
		% locally detrend the signal
		Y_ = Yi(id);
		% TODO, choose order
		if (o>0)
			A = vander_1d(cvec(id),o);
			b = A\Yi(id);
			Y_ = Y_ - A*b;
		end
		fyi = abs(fft(w.*Y_,ne)).^2;
		fy  = fy + fyi;
%		fy2 = fy2 + fyi.^2;
		fyi_(idx,:) = fyi;
	end
	% weighing factors, due to overlap at the window end points
	scale = 1.5/1.25; % TODO is this true?
	fy    = sqrt(fy);
	fyi_  = sqrt(fyi_);
	fys   = scale*std(fyi_)';
	% jackknife variance by leaving out one sample a time
%	fyi_ = sqrt(bsxfun(@minus,sum(fyi_),fyi_)*np/(np-1));
%	fys2 = (np-1)/np*bsxfun(@minus,fyi_',fy).^2;
%	fys  = sqrt(fys2);

	[fx Tx mask] = fourier_axis(Xi(1:ne));
	Tx  = Tx(mask);
	fy  = fy(mask);
	fys = fys(mask);
	% scale to true amplitude
	l   = length(X);
	fys  = fys/l;
	fy  = fy/l;
	fx  = fx(mask);
end


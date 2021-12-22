% 2016-02-12 14:51:25.412087030 +0100
% Karl Kastner, Berlin
%% autocorrelation function
% function acf = autocorr(y,numLags,correct_bias)
function acf = autocorr(y,numLags,correct_bias)
	if (isvector(y))
		y = cvec(y);
	end
	if (nargin<2 || isempty(numLags))
		numLags = size(y,1)-1;
	end

	y = bsxfun(@minus,y,mean(y));

	n = length(y);
	if (0)
		nFFT = 2^(nextpow2(n)+1);
	else
		nFFT = n;
	end
	%F = fft(y-mean(y),nFFT);
	F   = fft(y,nFFT);
	F   = F.*conj(F);
	acf = ifft(F);
	acf = acf(1:(numLags+1),:); % Retain non-negative lags
	%acf = acf./acf(1); % Normalize
	acf = real(acf);
	%acf = bsxfun(@times,acf,1./acf(1,:));
	acf = acf./acf(1,:);

	% correct for bias
	if (nFFT > n)
%nargin() > 2 && correct_bias)
		n = length(y);
		l = numLags+1;
		
		%size(y,1);
		acf = acf.*(n./(n-(1:l)'));
	end
end

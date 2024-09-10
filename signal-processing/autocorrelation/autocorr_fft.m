% 2016-02-12 14:51:25.412087030 +0100
% Karl Kastner, Berlin
%
%% estimate sample autocorrelation function
%
% function acf = autocorr(y,numLags,correct_bias)
function acf = autocorr_fft(y,numLags,correct_bias)
	if (isvector(y))
		y = cvec(y);
	end
	if (nargin<2 || isempty(numLags))
		numLags = size(y,1)-1;
	end

	y = bsxfun(@minus,y,mean(y));

	n = size(y,1);
	if (0)
		nFFT = 2^(nextpow2(n)+1);
	else
		nFFT = n;
	end
	%F = fft(y-mean(y),nFFT);
	F   = fft(y,nFFT);
	S   = F.*conj(F);
	acf = ifft(S);
	acf = acf(1:(numLags+1),:);
	acf = real(acf);
	% normalize
	acf = acf./acf(1,:);

	% correct for bias
	if (nFFT > n)
	%if (nargin() > 2 && correct_bias)
		n = length(y);
		l = numLags+1;
		acf = acf.*(n./(n-(1:l)'));
	end
end

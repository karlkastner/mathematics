% 2021-07-21 21:02:22.740297299 +0200
% Karl Kästner, Berlin
%
%% compute the normalized periodogram
%
% function S = periodogram(y,L)
function S = periodogram(y,L,nf)
	if (nargin()<2)
		L = 1;
	end
	if (isvector(y))
		y = cvec(y);
	end
	n = size(y,1);
	if (nargin()<3)
		nf = n;
	end
	energy = sum(y.*y);
	f     = fft(y,nf);
	if (nf < n)
		% truncation
		S  = 2*L./(nf*energy).*(f.*conj(f));
	else
		% padding with zeros
		S  = 2*L./(n*energy).*(f.*conj(f));
	end
	if (0 == energy(1))
		S = zeros(nf,1);
	end
end % periodogram


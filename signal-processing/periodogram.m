% 2021-07-21 21:02:22.740297299 +0200
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
	%var_y = var(y);
	var_y = rms(y).^2;
	%y  = y-mean(y);
	f     = fft(y,nf);
	if (nf < n)
		% truncation
		S  = 2*L./(n*nf*var_y).*(f.*conj(f));
	else
		% padding with zeros
		S  = 2*L./(n*n*var_y).*(f.*conj(f));
	end
	if (0 == var_y(1))
		S = zeros(nf,1);
	end
	%S  = f.*conj(f);
	%df = 1/L;
	%S  = 2*S./(df*sum(S));
end % periodogram


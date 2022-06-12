% 2021-07-21 21:02:22.740297299 +0200
% Karl Kästner, Berlin
%
%% compute the normalized periodogram in two dimensions
%
% function S = periodogram(y,L)
function S = periodogram_2d(y,L)
	if (nargin()<2)
		L = [1,1];
	end
	n      = size(y);
	energy = sum(sum(y.*y));
	f      = fft2(y);
	% padding with zeros
	S  = 2*L(1)*L(2)./(n(1)*n(2)*energy).*(f.*conj(f));
	if (0 == energy)
		S = zeros(n);
	end
end % periodogram_2d


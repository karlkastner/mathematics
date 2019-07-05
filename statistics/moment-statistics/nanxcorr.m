% Tue 12 Jul 17:08:47 CEST 2016
%
% cross correlation when values are nan
function x = nanxcorr(a,b,m)
	if (nargin() < 3)
		m = length(a);
	end
	x = zeros(2*m-1,1);

if (0)
	l = length(x);
	for idx=0:m-1;
		x(idx+m,1) = a(idx+1:end)'*b(1:end-idx)*(l/(l-idx));
		x(m-idx,1) = a(1:end-idx)'*b(idx+1:end)*(l/(l-idx));
	end
	x = 1/(sqrt(a'*a*(b'*b)))*x;
end
	l = length(a);
	for idx=0:m-1;
		% normalise with (l-0) to obtain biased estimate
		x(idx+m,1) = nansum(a(idx+1:end).*b(1:end-idx))/(l-idx);
		x(m-idx,1) = nansum(a(1:end-idx).*b(idx+1:end))/(l-idx);
	end
	% make first coefficient 1
	% note that the matlab xcorr 'coeff' estimate is biased
	% (closer to zero than actual cross correlation with increasing lags)
	x = 1/(nansum(a.*b/l))*x;

%for idx=-n:0
%	x(idx) = nancorr(a(idx:end-n+idx),b(n:end));
%end
%for idx=1:n
%	x(idx) = nancorr(a(1:end-idx),
%end
end

% 2017-02-16 14:40:30.783086307 +0100
%% filter with fast fourier transform
function [m s] = filtfft(f,p,fast)
	if (isvector(f))
		f=cvec(f);
	end
	n = size(f,1);
	id = (0:n-1)';
	l  = round(max(1,min(n/2,round(1/p*id))));
	r  = round(max(1,min(n/2,p*id)));
%	l  = round(max(1,min(n/2,round((1-p)*id))));
%	r  = round(max(1,min(n/2,(1+p)*id)));

	if (nargin()<3 || ~fast)
		m = NaN(size(f));
		s = NaN(size(f));

		for idx=1:n/2
			%l = max(1,round((1-p)*idx));
			%r = round(min(n/2,(1+p)*idx));
			m(idx) = mean(f(l(idx):r(idx),:));
			s(idx) = serr(f(l(idx):r(idx),:));
		end % for idx
	else
		% fast
		sf = cumsum(f);
		m  = (sf(r)-sf(l))./(r-l);
		s = [];
	end
end


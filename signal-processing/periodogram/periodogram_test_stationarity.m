% Wed  6 Oct 14:01:34 CEST 2021
%
%% test a periodogram for statoinarity
%% note : the method works, but is of little practical use,
%% as it requires about 50 periods and a small dx to detect a frequency change by a factor of 2
%
function [p,D,pp,ratio,SS,rr,mdx] = periodogram_test_stationarity(y,m,L,fmin,fmax)
%	y = y-mean(y);
	exclusive = true;
%	exclusive = false;
	y  = cvec(y);
	n  = size(y,1);
	ni = floor(n/m);
	if (nargin()<3)
		L = 1;
	end
	if (nargin()<3)
		fmin = 0;
	end
	if (nargin()<4)
		fmax = inf;
	end
	fx  = fourier_axis(L/m,ni);
	fdx = fx>=fmin & fx<=fmax; 
	nt = sum(fdx);
	for idx=1:m
		yi = y((idx-1)*ni+1:idx*ni);
		fi = fft(yi);
		SS(:,idx) = abs(fi).^2;
	end
if (exclusive)
	% mean excluding the respective slices
	S = 1/(m-1)*(sum(SS,2) - SS);
else
	S = mean(SS,2);
end
	%S = abs(fft(y)).^2;                                                             
	%S = S(1:end/2);                                                                 
	%SS = reshape(S,m,n/2/m)';
	%S = S(1:end/2);                                                                 
	%SS_ = reshape(S,m,n/m)';
	%S = m/n*mean(SS_,2);
	rr = SS(fdx,:)./S(fdx,:);
	ratio = rr;
	ratio = sort(flat(ratio));
	pp = (1:length(ratio))'/(length(ratio)+1);
if (exclusive)
	d1 = 2;
	d2 = 2*(m-1);
	pp(:,2) = fcdf(ratio,d1,d2);
	ratio(:,2) = finv(pp(:,1),d1,d2);
else
	a = 1;
	b = (m-1);
	pp(:,2) = betacdf(ratio/m,a,b);
	ratio(:,2) = m*betainv(pp(:,1),a,b);
end
	%D = max(abs(diff(pp,[],2)));
	[D,mdx] = max(abs(pp(:,2)-pp(:,1)));
	% p-value for rejecting the null hypothesis that the pattern is stationary
	% i.e. that the ratio does not follow the beta distribution
	p = 1 - kolmcdf(sqrt(nt).*D);
end


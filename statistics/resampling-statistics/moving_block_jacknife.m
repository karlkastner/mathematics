% Thu Apr 28 09:37:07 CEST 2016
% Karl Kastner, Berlin
%%
%% blocked Jacknfife for autocorrelated data
%% sliding block, statistically more efficient but computationally expensive
%% note, number of blocks must be sufficiently large h ~ sqrt(n)? << n
function s2 = moving_block_jackknife(fun,x,h)
	n = size(x,1);
	m = size(x,2);
	% get estimate for complete set
	Tn = fun(x);
	Tni = zeros(n-h,m);
	% get estimates with deleted blocks
	for idx=1:n-h+1
		% mask
		mask = true(n,1);
		mask(idx:idx+h-1) = false;
		% get estimate for current subset
		Tni(idx,:) = fun(x(mask,:));
	end
	% Dependence in Probability and Statistics, Bertail
%	s2 = h/((n-h)/(1)+1)*sum( (Tni - Tn).^2 );
	% bootstrap and jn, shao,1995, eq. 9.6, note, this misses bias correction for finite number of blocks / sample size
	% pseudo variable
	Tni_ = 1/h*bsxfun(@minus, n*Tn, (n-h)*Tni);
	% jacknife estimate of the (error) variance of the mean
	%s2 = h/(n*(n-h+1)) * sum( (Tni_ - 1/(n-h+1)*sum(Tni_)).^2 );
	s2 = h/(n*(n-h+1)) * sum( bsxfun(@minus,Tni_,1/(n-h+1)*sum(Tni_,1)).^2, 1);
	% unbiased variance estimate (this is 7.46 in Lahiri Resampling Methods for Dependent Data)
	s2 = s2*n/(n-h);
%	s2 = h/(n*(n-h+1)) * sum( bsxfun(@minus,Tni_,1/(n-h+1)*sum(Tni_,1)).^2, 1);

%	w = 1/(n-h+1);
%	w2 = 1/(n-h+1);
%	s2 = n/((n-h+1)) * w2 * sum( (Tni_ - 1/(n-h+1)*sum(Tni_)).^2 );

	% Encyclopedia of Environmetrics, Volume 1, el-shaarawi, 
	% pseudo variable
%	Tni_ = n*Tn - (n-h)*Tni;
	% jacknife estimate of the (error) variance of the mean
	%s2 = 1/(h*n*(n-h+1)) * sum( (Tni_ - 1/(n-h+1)*sum(Tni_)).^2 );
%	s2 = 1/(h*n*(n-h+1)) * sum( (Tni_ - 1/(n-h+1)*sum(Tni_)).^2 );
%	s2 = n*(n)/h^2*h/(n*(n-h+1)) * sum( (Tni - 1/(n-h+1)*sum(Tni)).^2 );
end


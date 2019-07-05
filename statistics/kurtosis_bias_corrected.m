% 2018-02-01 17:56:39.538255485 +0100
%
%% bias corrected kurtosis
function k = kurtosis_bias_corrected(x,resflag)
	if (isvector(x))
		x = cvec(x);
	end
	mu = mean(x);
	n  = size(x,1);
	if (nargin()<2 || ~resflag)
		e  = bsxfun(@minus,x,mu);
	else
		e = x;
	end
	% TODO this is not a good unbiased estimator of k2
	k2 = sum(e.^2)/(n-1);
	k = (n+1)*n/((n-1)*(n-2)*(n-3))*sum(e.^4)./k2.^2 - 3*(n-1).^2/((n-2)*(n-3)) + 3;
end


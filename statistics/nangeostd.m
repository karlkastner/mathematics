% 2017-06-28 13:03:02.637348502 +0200
%
%% geometric standard deviation ignoring nan-values
function x = nangeostd(x,dim);
	if (nargin()<2)
		dim = 1;
		if (isvector(x))
			x = cvec(x);
		end
	end
	mu = nangeomean(x,dim);
	x = exp(mean(log(bsxfun(@times,x,1./mu).^2)));
end


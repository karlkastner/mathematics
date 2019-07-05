% 2014-07-22 12:24:16.492622709 +0200
% Karl Kastner, Berlin
%
%% standard error of the mean of a set of uncorrelated samples
%
function s = serr(x,dim,cflag)
	if (isempty(x))
		s = NaN(class(x));
	elseif (nargin > 1 && ~isempty(dim))
		s = std(x,[],dim)./sqrt(size(x,dim));

	else
		if (isvector(x))
			s = sqrt(var(x)/length(x));
		else
			s = std(x,[],1)./sqrt(size(x,1));
		end
	end
		if (nargin() > 2 & cflag)
			res = x-mean(x);
			res = cvec(res);
			rho = res(1:end-1)'*res(2:end)/(res(1:end-1)'*res(1:end-1));
			f = f_correlation(rho,length(res));
			s = s*f;
		end
end


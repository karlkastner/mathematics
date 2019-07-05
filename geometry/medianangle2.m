% Sat 18 Feb 23:27:52 CET 2017
% Karl Kastner, Berlin
%
%% median angle
%%
%% input
%% alpha : x*m, [rad] angle
%%
%% ouput
%% ma    : 1*m, [rad] median angle
%% sa    : 1*m, [rad] standard error of median angle for uncorrelated error
%
% note: this is a multivariate median
function [ma, sa] = medianangle(alpha,dim)
	if (nargin()<2)
		dim = 1;
		if (isvector(alpha))
			alpha = cvec(alpha);
		end
	end
	n = size(alpha,1);

	s = sin(alpha);
	c = cos(alpha);

	ms = nanmedian(s,dim);
	mc = nanmedian(c,dim);

	% median angle
	ma = atan2(ms,mc);
	
	% residual
	if (nargout() > 1)
		res = bsxfun(@minus,alpha,ma);
		res = wrapToPi(res);

		% standard error
		% sa = serr(res);
		sa = sqrt(sum(res.^2)/(n*(n-1)));
	end
end % medianangle


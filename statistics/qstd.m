% 2015-04-06 15:32:33.996793669 +0200
% Karl Kastner, Berlin
% 
% Note : this is a measurement of dispersion and yields the same value for the
%        normal distribution as the standard deviation "std"
%        However, this is an own statistic and hence requires different
%        methods for calculating P-values and hypothesis testing
% TODO : is it more suitable to choose the quantiles at the inflection points
%	 of the normal dist or splitting the set in equal parts?
function sd = qstd(x,scaleflag,p) %,dim)
	if (nargin() < 3)
		p = 0.25;
	end
	p = min(p,1-p);
%	if (nargin() < 2 || scaleflag)
%		scale = 1./(2*norminv(1-p));
%	else
%		scale = 1;
%	end
	q = quantile(x,[p,1-p]);

	sd = qstdq(q,p);
end % function

%function s=qstd_(x)
%	c = normcdf(1);
%	q = quantile(x,[1-c c]);
%	if (isvector(x))
%		s = 0.5*(q(2)-q(1%));
%	else
%		s = 0.5*(q(2,:)-q(1,:));
%	end
%end

%function sd = qstd(x)
%	q = quantile(x,normcdf([-1 1]));
%	sd = 0.5*(q(2)-q(1));
%end


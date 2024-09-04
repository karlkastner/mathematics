% 2015-03-09 10:05:33.572003627 +0100
% Karl Kastner, Berlin
%
%% median and confidence intervals
%% c is a P value for the confidence interval,
%% default is 0.95 (2-sigma)
%% median of the colums of X
function [m, s, l, u] = median_man(X,P)
	if (isempty(X))
		m = NaN(class(X));
		s = NaN(class(X));
		l = NaN(class(X));
		u = NaN(class(X));
		return;
	end
	if (isvector(X))
		X = cvec(X);
	end
	if (nargin < 2)
		% P value
		% c = normcdf(1);
		%P = 0.95;
		% 1-sigma interval (~0.67)
			%P = (1-2*normcdf(-0.5))
		P = normcdf(1);
	end
	for idx=1:size(X,2)
%	n = size(X,2);
	n = sum(isfinite(X(:,idx)));
%	z  = norminv(0.5*c);
%	z = norminv(c);
	[pl pm pu] = median_ci(n,P);
%	qx    = quantile(X,[pl pm pu]);
%	qy    = quantile(Y,[pl pm pu]);
	q(:,idx) = quantile(X(:,idx),[pl pm pu]);
	end
	l = q(1,:);
	m = q(2,:);
	u = q(3,:);
	s = 0.5*(u-l);
end


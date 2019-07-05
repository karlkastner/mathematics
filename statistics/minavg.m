% 2014-09-01 21:14:18.277651012 +0200 minavg.m
%
%% solution of the minimum variance problem
%% minimise the variance of the weighted sum of n-independent
%% random variables with equal mean and individual variance
% for equal variances this is the average
% the solution is the harmonic mean (if variances are equal)
%
function [Y, Sy] = minavg(X,S)
	n = size(S,2);
	p = prod(S,2);
	den = zeros(size(S,1),1);
	q = zeros(size(S));
	for idx=1:n
		q(:,idx) = p./S(:,idx);
		den = den + q(:,idx);
	end
	coeff = zeros(size(S));
	for idx=1:n-1
		coeff(:,idx) = q(:,idx)./den;
	end
	coeff(:,end) = 1 - sum(coeff(:,1:end-1),2);

	Y =  sum(coeff.*X,2);
	if (nargout() > 1)
		Sy = sum(coeff.*S,2);
	end
%	Y =  coeff'*X;
%	if (nargout() > 1)
%		Sy = coeff'*S;
%	end
	%

	% denomminator = repmat(prod(S,2),1,size(S,2)-1);
	% nominator    = 

	% the optimal mean
	
	% the variance
end


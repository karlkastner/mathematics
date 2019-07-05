% Sat Nov 15 16:11:05 CET 2014
% Karl Kastner, Berlin
%
%% autocorrelation of a vector with nan-values
%
function [r rho] = autocorr1(X,n)
	if (nargin() < 2 | isempty(n))
		n = 1;
	end
	X = X(:);
	r = zeros(n,1);
	% problematic when it crosses nan values that separate columns
	for idx=1:n
		l = [NaN(idx,1); X(1:end-idx)];
		fdx = isfinite(l+X);
		r(idx,1) = (l(fdx)'*X(fdx))/(X(fdx)'*X(fdx));
	end
	for idx=1:n
		A(:,idx) = [NaN(idx,1); X(1:end-idx)];
	end
	fdx = isfinite(sum([A X],2));
	rho = A(fdx,:) \ X(fdx);
end


% Sat Nov 15 16:12:02 CET 2014
% Karl Kastner, Berlin
%
%% autocorrelation of a vector with nan-values
%
function [r roh] = autocorr2(X,n)
	if (nargin() < 2 | isempty(n))
		n = 1;
	end
	% TODO should be separate min for each dim
	n = min(n,min(size(X))-1);
%	left = [NaN(size(X,1),1)  X(:,1:end-1)];
%	top  = [NaN(1,size(X,2)); X(1:end-1,:)];
%	A = [top(:) left(:)];
%	fdx = isfinite(prod([A X(:)],2));
%	par = A(fdx,:) \ X(fdx);
	r = zeros(n,2);
	for idx=1:n
		l = [NaN(size(X,1),idx)  X(:,1:end-idx)];
		v = isfinite(l+X);
		r(idx,1) = l(v)'*X(v)/(X(v)'*X(v));
		t  = [NaN(idx,size(X,2)); X(1:end-idx,:)];
		v = isfinite(t+X);
		r(idx,2) = t(v)'*X(v)/(X(v)'*X(v));
	end
	if (nargout() > 1)
	for idx=1:n
		l = [NaN(size(X,1),idx)  X(:,1:end-idx)];
		left(:,idx) = l(:);
		t  = [NaN(idx,size(X,2)); X(1:end-idx,:)];
		top(:,idx) = t(:);
	end
	A = [top left];
	fdx = isfinite(prod([A X(:)],2));
	roh = A(fdx,:) \ X(fdx);
	roh = [roh(1:n) roh(n+1:end)];
	end
end


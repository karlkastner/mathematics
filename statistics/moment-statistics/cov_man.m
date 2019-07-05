% So 2. Aug 14:59:27 CEST 2015
% Karl Kastner, Berlin
%
%% covariance matrix of two vectors
%
function C = cov_man(x,y,w)
	mux = mean(x);
	if (nargin > 1 && ~isempty(y))
		muy = mean(y);
		dx  = x-mux;
		dy  = x-muy;
		if (nargin() < 3)
			n = length(dx);
			C = sum(dx.*dy)/(n-1);
		else
			C = sum(w.*dx.*dy)/(sum(w)-1);
		end
	else
		if (nargin() < 3)
			dx  = bsxfun(@minus,x,mux);
			n = length(dx);
			C = dx'*dx/(n-1);
		else
			mux = wmean(x,w);
			dx  = bsxfun(@minus,x,mux);
			W   = diag(sparse(w));
			C   = dx'*W*dx/(sum(w)-1);
		end
	end
end % cov_man


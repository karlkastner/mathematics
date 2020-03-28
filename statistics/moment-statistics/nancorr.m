% 2014-10-28 15:28:45.414260739 +0100
% Karl Kastner, Berlin
%
%% (co)-correlation matrix when samples a NaN
function [R, p] = nancorr(A,varargin)
	C = nancov(A,varargin{:});
%	if (nargin() < 2)
%		B = A;
%		C = nancov(A);
%	else
%		C = nancov(A,B,varargin{:});
%	end
	if (all(isfinite(C(:))))
		R = corrcov(C);
	else
		R = NaN(size(C));
	end
%corrcoef(a,b,'rows','pairwise')

%	fdx = isfinite(A.*B);
%	R = corr(A(fdx),B(fdx));
end

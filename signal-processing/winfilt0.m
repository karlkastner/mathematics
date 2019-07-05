% Mo 1. Sep 18:27:30 CEST 2014
% Karl Kastner, Berlin
%% filter with window
function [y err] = winfilt(x,lf,p,recursive,wtype,lower)
	if (nargin() < 4 || ~recursive)
		rmax = 1;
	else
		rmax = 10;	
	end
	if (nargin() < 5)
		wtype = 0;
	end
	if (nargin() < 6)
		lower = 0.5;
	end

	% allocate memory
	y   = zeros(size(x));
	err = zeros(size(x));
	for idx=1:length(x)
		l = max(1,idx-lf);
		r = min(idx+lf,length(x));
		% set up local regression matrix
		x_ = (l:r)'-idx;
		b  = x(l:r);
		% find nan input
		fdx = find(isfinite(b));
		b = b(fdx);
		if (isempty(fdx) | length(fdx)/(2*lf+1) < lower )
			y(idx) = NaN;
			err(idx) = NaN;
			continue;
		end
		A = [];
		for jdx=1:p+1
			A(:,jdx) = x_(fdx).^(jdx-1);
		end % for jdx

		% filter window
		if (wtype == 1)
		 	f  = (1-sin(0.5*pi*x_/(lf+1)).^2);
		else
			f = ones(size(x_));
		end
		f=f(fdx)/sum(f(fdx));

		W = diag(f);

		res = ones(size(f));
		for rdx=1:rmax
		% build new weight according to the residual (TODO, zero residual may lead to divergence)
		r2m = nanmedian(res.^2);
		WW = diag(1./max(r2m,res.^2))*W;
		% solve for coefficient
		c = (A'*WW*A) \ (A'*WW*b);
		% evaluate polynomial
		y(idx) = c(1);
		% residual
		res = A*c - b;

		end

		% residual norm
		% kish's weight estimate could be used, but here we assume half band width sample length
		%err(idx) = 1/sqrt(0.5*length(fdx))*sqrt(res'*res);
		err(idx) = sqrt(res'*W*res);
	end % for idx
end % winfilt


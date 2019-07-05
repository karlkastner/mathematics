% Mon Dec  1 14:22:20 CET 2014
% Karl Kastner, Berlin
%% short time fourier transform
% TODO use a window
%
function [obj] = transform(obj,time,val) 

	val = double(val);
	if (isvector(val))
		val = cvec(val);
	end


	% set up regression matrix
	[A1, A2, Dl, Dr] = obj.stftmat(time);
	obj.A1 = A1;
	obj.A2 = A2;

	% remove invalid rhs
	fdx1 = isfinite(val);
	fdx2 = flat(repmat(rvec(fdx1),obj.ni,1));

	time = time(fdx1);
	val  = val(fdx1);

	A1_  = A1(fdx1,fdx2);
	A2_  = A2(fdx2,:);

	% local regularisation, because:
	% the regression matrix will be ill conditioned
	% - at the end of the time series
	% - at gaps in the time series
	% - if the interval is too small
	% - if the 1./T exceeds sample frequencies

	A    = A1_*A2_;
	AA   = A'*A;
	d    = diag(AA);
	o    = sum(abs(AA),2)-d;
	fdx_ = (0 == d & 0 == o);
	w    = max(0,o-d);
	w(fdx_) = 1;
	DD   = Dl'*diag(sparse(w))*Dl + Dr'*diag(sparse(w))*Dr;

	% solve least squares problem
	if (0)
		fun     = @(c) A2_'*(A1_'*(A1_*(A2_*c)));
		c_      = minres(fun,(A2_'*(A1_'*val)),[],obj.maxit);
	else
		% about 6 times faster
		lambda = obj.lambda;
		%c_ = minres(A_'*A_,A_'*val_,[],obj.maxit);
		%c_ = minres(A_'*A_ + lambda*(D_'*D_),A_'*val_,[],obj.maxit);
		%condest(A_'*A_ + lambda*(D_'*D_))
		% TODO use chol
		c_ = (AA + DD) \ (A'*val);
	end

	c       = NaN(size(A2,2),1);
	c(~fdx_) = c_(~fdx_);
	c       = reshape(c,obj.ni,obj.nti);
	obj.coeff = c;



	if (0)

	nc = size(val,2);
	nt = size(val,1);

	% samples per interval
	n = round(Ti/dt);

	% number of intervals
	m = floor(nt/n);

	% time at interval centres
	tc     = n*dt*((0:m-1)'+0.5);
	timei  = dt*(0:n-1)';

	A = fouriermtx(timei,T);

	c = zeros(1+2*length(T),m,nc);
	for cdx=1:size(val,2)
		% split into intervals, discard tail
		vali = reshape(val(1:m*n,cdx),n,m);

		% regress constituent coefficients
		c(:,:,cdx)      = A \ vali;
	end %  cdx
%	res    = A*c - vali;
%	serr   = NaN;
	

	obj.tc    = tc+t0;
	obj.coeff = c;
	end

%	serr   = sqrt(diag(res.*res)/(n-2*length(T)-1));
end % stft

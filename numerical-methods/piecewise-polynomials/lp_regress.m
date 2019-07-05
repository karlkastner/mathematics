% Wed Mar  4 17:45:14 CET 2015
% Karl Kastner, Berlin
%
% fit pieceweise polynomials with lagrangian basis
%
function [c, serr, cnt, sdx] = lp_regress(order,x0,x,val)
	if (isempty(order))
		order = 1;
	end
	% TODO : so far order has to be constant or linear
	order = max(0,min(1,order));
	n   = length(x0);
	dof = n;
%	dof = n+1;
	cnt = zeros(n,1);
	c    = NaN(dof,1,class(val));

	% filter
	fdx = isfinite(val);
	x   = x(fdx);
	val = val(fdx);

%	% TODO this should be combined with the QR factorisation
%	% the current QR implementation fails for l(x) < dof
	if (length(x)<dof)
		serr = NaN(class(val));
		return;
	end


	%dx = (x0(2)-x0(1));
	dx = (x0(end)-x0(1))/(length(x0)-1);
%	dx = (x0(end)-x0(1))/(length(x0));
	% determine the leftmost support point for each source point
	%x_ = x - 0.5*dx*order;
	x_  = x;
	x0_ = x0 - 0.5*dx+0.5*dx*order;
	% TODO, if the grid is equally spaced, this is simply floor((t-t0)/dt)
	sdx = sum(bsxfun(@gt,cvec(x_),rvec(x0_)),2);


	% limit
	sdx = min(n,max(1,sdx));

	% determine rightmost index
	rdx=sdx+order;
	rdx=min(n,max(1,rdx));

	% ensure sufficient points at boundary
	rdx(sdx==1) = order+1;
	sdx(rdx==n) = n-order;

%figure()
%plot(x,[1+x*(n-1) sdx rdx])
%pause


	% count (TODO only valid for order==1)
	if (nargout() > 2)
		for idx=1:n
			cnt(idx) = sum(sdx == idx);
		end
	end

	% shift and scale source point coordinates to unit interval
	x = (cvec(x) - cvec(x0(sdx)))./dx; %cvec(x0(sdx+1)-x0(sdx));
%	x = x-order/2;
%	x = x-1/2;
%x

	% TODO, why limt?
	  % x = min(1,max(0,sdx));

%	% find and average identical points
%	[x id ai] = unique(x);
%	val = accumarray(ai,val) ./ accumarray(ai,ones(size(val)));

	% vandermonde matrix of source points coordinates
	Ax = vander_1d(x,order);

	% inverse of the Vandermonde matrix on the unit interval
	% the leftmost point is the origin
	A0i = inv(vander_1d((0:order)',order));
%	switch (order)
%	case {0}
%		A0i = 1;
%	case {1}
%		A0i = [ 1 0;
 %       	       -1 1];
%	case {2}
%		A0i = inv(vander_1d([0; 1; 2],2));
%	otherwise
%		error('not yet implemented')
%	end

	% transformation of sample to basis points
	C = Ax*A0i;

	% construct regression matrix
	% A sparse matrix could be used here, however, a full matrix is optimal
	% if the number of source points is larger than the number of segments,
	% which has to be the case to make the regression matrix non-singular
	A = zeros(length(val),dof);
	for idx=1:length(val)
		A(idx,sdx(idx):rdx(idx)) = C(idx,:);
	end
	% only compute values at segments at points that are supported,
	% i.e. reduce A to non-singular sub-matrix
	% rather than the pseudo inverse one of the linear dependend target points are excluded
	[Q R]    = qr(A,0);
	% TODO this is also bogus, use svd or set to nan
	valid    = abs(diag(R)).^2 > eps(class(val));
	opt.UT   = true;
	c(valid) = linsolve(R(valid,valid),Q(:,valid)'*val,opt);
	
	%valid = any(A)';
	% determine values at basis points
	% (polynomial coefficients are indicrectly determined)
	% c(valid) = A(:,valid) \ val;
	% standard error
	res = A*c - val;
	serr = sqrt(res'*res/(length(val)-length(c)));
end % lp_regress


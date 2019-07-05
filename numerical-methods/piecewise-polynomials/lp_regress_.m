% Wed Mar  4 17:45:14 CET 2015
% Karl Kastner, Berlin
%
% fit piecewise polynomials with lagrangian basis
%
% n : number of segments
function [c, serr] = lp_regress(order,x0,x1,n,x,val)
	% TODO : so far order has to be constant or linear
	order = max(0,min(1,order));
	n = max(1,n);
	switch(order)
		case {0}
			dof = n;
		case {1}
			dof = n+1;
		otherwise
			error('lp_regress');
	end

	% filter
	fdx = isfinite(val);
	x   = x(fdx);
	val = val(fdx);

%	% TODO this should be combined with the QR factorisation
%	% the current QR implementation fails for l(x) < dof
	if (length(x)<dof)
		c    = NaN(dof,1,class(val));
		serr = NaN(class(val));
		return;
	end

	% segment length
	dx = (x1-x0)/n;
	% determine segment for each source point
	sdx = floor((x-x0)/dx)+1;
	sdx = min(n,max(1,sdx));
	% shift and scale source point coordinates to unit interval
	x = (x - x0)/dx - sdx + 1;

%	% find and average identical points
%	[x id ai] = unique(x);
%	val = accumarray(ai,val) ./ accumarray(ai,ones(size(val)));

	% vandermonde matrix of source points
	Ax = vander_1d(x,order);
	% inverse of the Vandermonde matrix on the unit interval
	switch (order)
	case {0}
		A0i = 1;
	case {1}
		A0i = [ 1 0;
        	       -1 1];
	otherwise
		error('not yet implemented')
	end
	% transformation of sample to basis points
	C = Ax*A0i;
	% construct regression matrix
	% A sparse matrix could be used here, however, a full matrix is optimal
	% if the number of source points is larger than the number of segments,
	% which has to be the case to make the regression matrix non-singular
	A = zeros(length(val),dof);
	for idx=1:length(val)
		A(idx,sdx(idx):sdx(idx)+order) = C(idx,:);
	end
	% only compute values at segments at points that are supported,
	% i.e. reduce A to non-singular sub-matrix
	% rather than the pseudo inverse one of the linear dependend target points are excluded
	[Q R]    = qr(A,0);
	valid    = abs(diag(R)).^2 > eps(class(val));
	c        = NaN(dof,1,class(val));
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


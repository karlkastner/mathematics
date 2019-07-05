% Sun Feb  1 17:31:40 CET 2015
% Karl Kastner, Berlin
%
%% lagrangian basis piecwie interpolation, predicor
function val = lp_predict(order,x0,c,x)
	n = length(x0);
	c = cvec(c);

	% degrees of freedom (number of basis points)
	dof = length(c);

	dx = (x0(end)-x0(1))/(length(x0)-1);
%	dx = (x0(end)-x0(1))/(length(x0));
%	dx = (x0(2)-x0(1));

	% determine segment index of target points
	% determine the leftmost support point for each source point
	%x_ = x - 0.5*dx*order;
	x_ = x;
	x0_ = x0 - 0.5*dx+0.5*dx*order;
	% TODO, if the grid is equally spaced, this is simply floor((t-t0)/dt)
	sdx = sum(bsxfun(@gt,cvec(x_),rvec(x0_)),2);
	% sdx = sum(bsxfun(@gt,cvec(x),rvec(x0)),2);

	% limit
	sdx = min(n,max(1,sdx));

	% determine rightmost index
	rdx=sdx+order;
	rdx=min(n,max(1,rdx));

	% ensure sufficient points at boundary
	rdx(sdx==1) = order+1;
	sdx(rdx==n) = n-order;

	% dx = cvec(x0(sdx+1)-x0(sdx));
	% shift and scale target point coordinates to local unit interval
	x = (cvec(x)-cvec(x0(sdx)))./dx;
%	x = x-order/2;
%	x = x-1/2;

	% limit
%	x = min(1,max(0,x));

	% vandermonde matrix of source points coordinates (x)
	Ax = vander_1d(x,order);

	% inverse of the Vandermonde matrix on the unit interval
	% TODO the lagrange polynomial has a direct solution formula for the inverse
	A0i = inv(vander_1d((0:order)',order));
%	switch(order)
%	case {0}
%		A0i = 1;
%	case {1}
%	        A0i = [ 1, 0
 %       	       -1, 1];
%	case {2}
%		
%	otherwise
%		error();
%	end

	% transformation of support to targe point values
	C = Ax*A0i;

	% predict
	val = zeros(size(x));
	%A = zeros(length(val),dof);
	for idx=1:length(val)
		% implicite evaluate matrix vector product
		val(idx) = C(idx,:)*c(sdx(idx):sdx(idx)+order);
		%A(idx,sdx(idx):rdx(idx)) = C(idx,:);
	end
	% TODO prediction error
end % lp_predict


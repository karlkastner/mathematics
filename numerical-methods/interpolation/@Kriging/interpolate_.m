% Sun Jun 15 21:46:22 WIB 2014
% Karl Kastner, Berlin
%
%% interpolate with Krieging method
%%
%% this function may interpolate several quantities per coordinate,
%% using the same variogram, if the semivariance of the quantities differs,
%% the user may prefer to estimate the semivariance and interpolate each quantity
%% individually
%%
%% Xs  : source point coordinates
%% Vs  : value at source points
%% Xt  : targe point coordinates
%% Vt  : value at target points
%% E2t : squared interpolation error at target points
%
% TODO set radius according to semivariance parameter
function [vt, et, ot, jd, res, dx, dy, obj] = interpolate_(obj, Xs, Vs, Xt, Rmin, func)
	res = [];
	dx  = [];
	dy  = [];
	nt  = lenght(Xt);

	% shift points to avoid cancellation errors
	% TODO set Xt to zero in mind only
	for idx=1:nt
		Xs(:,idx) = Xs(:,idx) - Xt(idx);
		Xt(idx) = 0;
	end

	% distance from measured to interpolated point
	% euclidean, not mahalanobis distance
	for idx=1:nt
		Db(:,idx) = Xs(:,idx) - Xt(idx);
	end
	Db2 = Db(:,1).^2 + Db(:,2).^2;

	% select indices where the variance is smaller 0.95
	lt = log(1-obj.thresh);
	switch (obj.model.get())
		case {Model.EXPONENTIAL}
			% exponent of the semivariance
			par = obj.param(2:end);
			E = abs(Db)*par(:);
			% exp(-E) < 0.95 <=> -log(0.95) < E
			jd = find( -lt > E & Db2 > Rmin*Rmin );
		case {Model.GAUSSIAN}
			% exponent of the semivariance
			E = (Db.^2)*obj.param(2:end);
			% exp(-E) < 0.95 <=> -log(0.95) < E
			jd = find( -lt > E & Db2 > Rmin*Rmin);
		otherwise
			error();
	end % switch

	% invalidate out of range points
	if (isempty(jd))
%'molch'
%pause
		vt  = NaN;
		et = NaN;
		ot = NaN;
	else
	
	One = ones(length(jd),1);

	% distance for each spatial dimension
	% this is the euclidean, not the mahalanobis distance
	% scale coordinates system
	%Xs = Xs/obj.Rmax; -> only for polynomial
	for idx=1:nt
		Da(:,:,idx) = Xs(jd,idx)*One' - One*Xs(jd,idx)';
		Db(jd,idx)  = Xs(jd,idx) - Xt(idx);
	end % for idx

	% set up of semivariance Kriging matrix
	Dav = reshape(Da,[],nt);
	A_ = obj.model.svfunc(obj.param, Dav);
	A_ = reshape(A_,length(jd),length(jd));
	% set up the Kriging matrix
	order = obj.order.get();
	% minimum number of points for polynomial of order n is:
	% 1/fact(ndims) \prod_{i=1}^ndims (n+i)
	%xdims = length(Xt);


	if (Order.CONSTANT == order || length(jd) < dof(nt,1))
		A = [  A_  One;
		      One'   0];
		% set up regression vector
		b = [ obj.model.svfunc(obj.param, Db(jd,:)); 1 ];
		ot = 0;
	elseif (Order.LINEAR == order || length(jd) < dof(nt,2))
		% TODO maybe transform coordinates with Xt = 0, to avoid cancellation of digits
		% TODO, to estimate the coefficients, at least 3 input points are required
		% so better switch to constant interpolation for this case
		A = [ A_        , One, Xs(jd,:);
		      [One'; Xs(jd,:)'], zeros(nt+1) ];
		% set up regression vector
		b = [ obj.model.svfunc(obj.param, Db(jd,:)); 1; Xt(:) ];
		ot = 1;
	elseif (Order.QUADRATIC == order)
		% TODO this is only for 2d at the moment
		% set up vandermonde matrix
		vander = [One Xs(jd,:) Xs(jd,:).^2 Xs(jd,1).*Xs(jd,2)];
		A = [     A_,   vander
		     vander', zeros(6)];
		% set up regression vector
		b = [ obj.model.svfunc(obj.param, Db(jd,:)); 1; Xt(:); Xt(:).^2; Xt(1).*Xt(2) ];
		ot = 2;
	else
		error('');
	end % if

	% regress Lagrange multipliers
	if (length(b) < obj.MAX_SIZE_DIRECT)
		% use direct solver for small matrices
		w = A\b;
	else
		% use iterative solver
		% default number of iterations is just 10,
		% but roughly n = N^1/dim iterations are at least necessary for convergence
		%maxit = ceil(min(size(A))^(1/length(Xt)));
%		w = gmres(A,b,[],obj.TOL,maxit);
		maxit = length(jd);
		w = minres(A, b, obj.TOL, maxit);
		% TODO, check convergence
	end % if

	% estimate the value
	vt = w(1:length(jd))'*Vs(jd,:);

	% squared error of the estimate
	% this is also true for linear and quadratic interpolation
	% as the target point coordinates are in the origing (Xt = y0 = 0)
	% cf Computing Risk for Oil Prospects: Principles and Programs: Principles, p169 	
	% Atlas of Antarctica: Topographic Maps from Geostatistical Analysis of, 52
	% Error Propagation in Environmental Modelling with GIS, p21
	% s^2 = b'*A^-1*b 
	et = w'*b;
%'honk'
%pause

%	end % sufficient number of points in proximity
%	end % not idw
%	end % for each target point
end % interpolate


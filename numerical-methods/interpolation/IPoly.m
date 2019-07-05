% Sun Jul  6 16:59:24 WIB 2014
% Karl Kastner, Berlin
%% polynomial interpolation class
% TODO infomogeneous class design, design access to cfunc and svfunc (Kriging) in a similar manner
% TODO support anisotropic calculation of the distance
% TODO cancel source points, where error is large (can be done together with verification)
% kriging : clustering before, only n nearest neighbours into matrix (sparse matrix solution)
classdef IPoly < Interpolator
	properties
		%Rmax
		%Rmin
		%order
		%nverify

		% estimated covariance
		ecov
		% covariance / weighing function
		cfunc
		%rtol
		%s
	end % properties
	methods
	function obj = IPoly(Rmax,Rmin,order,nverify,aspect_ratio)
		% call superclass constructor
		obj = obj@Interpolator(Rmax,Rmin,order,nverify,aspect_ratio);
		obj.ecov.n  = 100; 
	%	obj.ecov.C  = @cfunc;
		obj.ecov.param  = [1 1];
		addpath([ROOTFOLDER,fielsep(),'../master/thesis/src/sandbox/ported/']); 
	end % constructor

	function obj = estimate_parameter(obj,Xs,Vs,Xt)
		% if this is done with the source points the minimum is zero!!!
		% the cut-off radius is somewhat benign, here a quick quadratic estimation
		% to minimise the verification error
		nverify = obj.nverify;
		% estimate to 5% accuracy
		obj.nverify = round(1/(0.05^2));
		Rmax = obj.Rmax*[0.25 0.5 1 1.5 2 2.5 3 ];
		%Rmax = obj.Rmax*[0.5 1 1.5];
		obj.Rmin = 0;
		% TODO should be nopt
		vdx = randi(size(Xt,1),obj.nverify,1);
		for idx=1:length(Rmax)
			% verify is no good with positive radius, as the error goes to zero for cutoff 0
			% so the interpolation error estimate 
			obj.Rmax = Rmax(idx);
			%obj.Rmin = 0.5*Rmax(idx)
			%[Es vdx] = verify(obj,Xs,Vs);
			%q = quantile(Es,[normcdf(-1),normcdf(+1)]);
			%E(idx) = 0.5*(q(2)-q(1));
			[val Es] = obj.interpolate(Xs,Vs,Xt(vdx,:),0);
			%E(idx) = nanmedian(Es);
		end % for idx
		% set up quadratic polynomial
		p3 = polyfit(Rmax,E,2);
		% derive
		p2 = [2*p3(1) p3(2)];
		% find zero
		Ropt = -p2(2)/p2(1);
		fprintf('Optimum cut off radius is %f\n',Ropt);
		% checks
		if (Ropt <= 0)
			error('Optimum out of bounds');
		end
		% derive another time
		p1 = p2(1);
		if (p1 <= 0)
			error('Function not convex around initial condition, no minimum found');
		end
		% set the optimal parameter
		obj.Rmax = Ropt;
		obj.Rmin = 0.5*Rmax;
		% reset nverify
		obj.nverify = nverify;
		plot(Rmax,E,'k.'); x=linspace(0,1000)'; hold on; plot(x,[x.^2 x x.^0]*p3','k');
	end % estimate_parameter

	function obj = estimate_covariance(obj,Xs,Vs,Xt)
		disp('Estimating covariance model parameters');
		n = obj.ecov.n;
		% select a subset of points
		id = randi(size(Xt,1),n,1);
%		cfunc = @(dx,dy) speye(length(dx));
%		cfunc = @(dx,dy,
		Res = [];
		Dx  = [];
		Dy  = [];
		for idx=1:n
			x0 = Xt(id(idx),:);
			% param should be 0, but this will lead to overflow
			%func = @(dx,dy) feval(@cfunc,dx,dy,x0(1),x0(2),obj.ecov.param);
			% func = @(dx,dy) eye(length(dx));
			%func = @(dx,dy
			[vt et ot fdx res dx dy] = obj.interpolate_(Xs,Vs,x0,0,func);
			Res = [Res;res];
			Dx  = [Dx; dx];
			Dy  = [Dy; dy];
		end % for idx
		% remove outliers
		fdx = find(Res.*Res < quantile(Res.*Res,0.95));
		Res = Res(fdx);
		Dx = Dx(fdx);
		Dy = Dy(fdx);
		% regress the covariance parameters
		s2 = Res'*Res/length(Res);
		b = log((1 - Res.*Res/s2).^2 );
		A = [-abs(Dx) -abs(Dy)];
		param = A \ b;
		obj.ecov.s2     = s2;
		obj.ecov.param  = param;
		obj.ecov.id     = id;
	%	x = linspace(0,750,10);
	%	B = bin2d(abs(Dx),abs(Dy),Res.^2,x,x);
	%	bar3(x,B)
	end % estimate covariance

	function [vt et ot fdx res dx dy obj] = interpolate_(...
						obj,Xs,Vs,x0,Rmin) %,cfunc)
		% find all points within radius Ri
		% TODO use a quadtree here
		% TODO is bsxfunc faster?
		Rmax2 = obj.Rmax*obj.Rmax;
	
%		if (1)
%			D = obj.dist(Xs,repmat(x0,size(Xs,1),1));
%			D2 = (D(:,1).^2 + D(:,2).^2);
%			fdx = find(D2 < Rmax2 & D2 >= Rmin*Rmin);
%			D = D(fdx,:);
%			fdx_ = fdx;
%		else
			gdx = obj.qtree.box_neighbour(x0(1),x0(2),max(1./obj.s)*obj.Rmax);
			D   = obj.dist(Xs(gdx,:),repmat(x0,size(gdx,1),1));
			D2  = (D(:,1).^2 + D(:,2).^2);
			fdx = find(D2 < Rmax2 & D2 >= Rmin*Rmin);
			D   = D(fdx,:);
			fdx = gdx(fdx);
%		end
%		if (length(fdx_) ~= length(fdx) || sum((sort(fdx(:))-sort(fdx_(:))).^2) > 0)
%			error()
%		end
		dx = D(:,1);
		dy = D(:,2);
		% TODO, check number for degree of poly
		order = obj.order;
		while (order >= 0)
			nmin = dof(size(Xs,2),order);
			if (length(fdx) >= nmin)
				break;
			end
			order=order-1;
		end %while
		if (order < 0)
			vt = NaN;
			et = NaN;
			ot = NaN;
			dx = [];
			dy = [];
			res = [];
			return;
		end
		% set up regression matrix
		% the coordinates are shifted with respect to x0,y0
		% and scaled to 1 to avoid cancellation error
		A = vander_2d(D/obj.Rmax,order);
		b = Vs(fdx,:);
		W = feval(obj.cfunc,dx/obj.Rmax,dy/obj.Rmax);
		% weighted least squares matrix
		%wii = obj.Rmax-sqrt(D2(fdx));
		%wii = wii/sum(wii);
		%W = spdiags(wii,0,length(fdx),length(fdx));
		% regress local coefficients
		Wi = inv(W);
		while (rcond(full(A'*Wi*A)) < obj.rtol)
			% degcrease order
			order = order-1;
			if (order < 0)
				vt = NaN;
				et = NaN;
				ot = NaN;
				res = [];
				dx  = [];
				dy = [];
				return;
			end %if
			% lazy, existing matrix could also be cut	
			A = vander_2d(D/obj.Rmax,order);
		end  %while
		coeff = (A'*(Wi*A)) \ (A'*(Wi*b));
		% evaluate at target point, only constant term is nonzero, as point in origin
		vt = coeff(1);
		% estimation error
		%wii = spdiags(W,0,length(W),lengt(W));
		%wii  = wii/sum(wii);
		%Et(idx,1) = norm(A*coeff-b)/sqrt(length(fdx)-length(coeff));
		res  = (A*coeff-b);
		n    = length(res);
		ndof = length(coeff);
		et   = sqrt(n/(n-ndof)*(res'*(Wi*res)));
		ot = order;
	end % interpolate_
	end % methods
end % classdef


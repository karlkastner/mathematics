% Sun Jan 11 13:30:06 CET 2015
% Karl Kastner, Berlin
%% sparse interpolation class
% TODO make this function compatible for 1D,2D,3D at the moment only 2d
classdef ISparse < Interpolator
	properties
		% number of points
		np
	end
	methods
	function obj = ISparse(Rmax,Rmin,order,nverify,aspect_ratio,np);
		% call superclass constructor
		obj = obj@Interpolator(Rmax,Rmin,order,nverify,aspect_ratio);
		% set degree of freedom
		mindof = (obj.order+1).*(obj.order+2)/2 + 1;
		if (mindof > np)
			str = strpintf('Increased number of points to %d, the minimum number of points for $d-order interpolation and error estimation');
			warning(str);
			obj.np = mindof;
		else
			obj.np = np;
		end
	end

	function [vt et ot id res dx dy obj] = interpolate_(obj,Xs,Vs,x0,Rmin,cfunc)
	id = zeros(obj.np,1);
	% distance from all source to current target point
	dx = Xs(:,1) - x0(1);
	dy = Xs(:,2) - x0(2);
	switch (obj.order)
	case {0} % constant interpolation, minimise linear term
		id = quartant_selection(dx,dy);
%{
		% this is a greedy point selection algorithm
		for idx=1:obj.np
			% residual weight of selected points
			rdx = dx + sum(dx(id(1:idx-1)));
			rdy = dy + sum(dy(id(1:idx-1)));
			% norm of residual distances
			% by independce square of sum is sum of squares
			% this is optimal for isotropic functions
			res = rdx.^2 + rdy.^2;
			% select point that minimises the error
			[mv id(idx)] = min(res);
		end
%}
		% interpolation value
		vt = mean(Vs(id));
		% residual
		res = vt-Vs(id);
		% standard error
		et = sqrt(res'*res/(length(id)-1));
		%et = sqrt(res'*res/(obj.np-1));
		ot = obj.order;
	case {1}
		%id_ = quartant_selection(dx,dy);
		id = obj.qtree.nearest_quad(x0(1),x0(2));
		%if (length(id) ~= length(id_) || norm(sort(id)-sort(id_)) > 0)
		%	%error
		%	'warning'
		%	id
		%	id_
		%end

		% linear interpolation, minise sum of quadratic and cubic term
		% cubic term needs to be includes as quadratic term is unsigned
%		dx2 = dx.*dx;
%		dy2 = dy.*dy;
%		dx3 = dx2.*dx;
%		dy3 = dy2.*dy;
%		flag = zeros(size(dx));
%		for idx=1:obj.np
			% residual weight of selected points
%			rdx2 = dx2 + sum(dx2(id(1:idx-1)));
%			rdx3 = dx3 + sum(dx3(id(1:idx-1)));
%			rdy2 = dy2 + sum(dy2(id(1:idx-1)));
%			rdy3 = dy3 + sum(dy3(id(1:idx-1)));
%			res = rdx2.^2 + (1/3*rdx3).^2 + rdy2.^2 + (1/3*rdy3).^2;
			%res2 = (0.5*(dx2 + sum(dx2(id(1:idx-1))))).^2 ...
			%+ (1/6*(dx3 + sum(dx3(id(1:idx-1))))).^2 ...
			%+ (0.5*(dy2 + sum(dy2(id(1:idx-1))))).^2 ...
			%+ (1/6*(dy3 + sum(dy3(id(1:idx-1))))).^2;
%			res2 = ... %(0.5*(dx2 + sum(dx2(id(1:idx-1))))).^2 ...
%			+ (1/6*(dx3 + sum(dx3(id(1:idx-1))))).^2 ...
%			... %+ (0.5*(dy2 + sum(dy2(id(1:idx-1))))).^2 ...
%			+ (1/6*(dy3 + sum(dy3(id(1:idx-1))))).^2 ...
%			+ flag;
%			[mv id(idx)] = min(res2);
%			flag(id(idx)) = NaN;
%		end
		if (length(id) < 1)
			vt = NaN;
			et = NaN;
			ot = NaN;
			return;
		end

		if (length(id) > 2)
			% 2d vandermonde matrix, x0 shifted to origin
			A = [ones(length(id),1) dx(id) dy(id)];
		else
			A = ones(length(id),1);
		end
		% interpolation coefficient
		c = A \ Vs(id);
		% interpolation value (polynomial evaluated at origin)
		vt = c(1);
		% residual
		res = A*c - Vs(id);
		% standard error
		et = sqrt(res'*res/(obj.np-3));
		otherwise
			error('Not yet implemented')
	end % case
	% The order does not vary
	ot = obj.order;
	fdx = id;
	end % interpolate_
	end % methods
end % classdef


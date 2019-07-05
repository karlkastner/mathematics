% Sun Jul  6 12:02:52 WIB 2014
% Karl Kastner, Berlin
%
%% nearest neighbour interpolation
classdef Inn < Interpolator
	properties
		% no additional properties
	end % properties

	methods
	function obj = Inn(Rmax,Rmin,order,nverify,aspect_ratio)
		% call superclass constructor
		obj = obj@Interpolator(Rmax,Rmin,order,nverify,aspect_ratio);
	end

	function [vt et ot fdx res dx dy obj] = interpolate_(obj,Xs,Vs,x0,Rmin,cfunc)
		% find all points within radius Ri
		% TODO use a quadtree here
		% TODO is bsxfunc faster?
		D = obj.dist(Xs,repmat(x0,size(Xs,1),1));
		D2 = (D(:,1).^2 + D(:,2).^2);
		R2max = obj.Rmax*obj.Rmax;
		R2min = Rmin*Rmin;
		[d2min fdx] = min(D2+realmax('double')*(D2<R2min));
		%dx = D(fdx,1);
		%dy = D(fdx,2);
		if (d2min < R2max)
			% choose value of nearest neighbour
			vt = Vs(fdx);
			% TODO, error estimate
			et = NaN;
			% actually the order is not well defined
			ot = 0;
		else
			vt = NaN;
			et = NaN;
			ot = NaN;
		end
		res = [];
		dx = [];
		dy = [];	
		% constant interpolaton
		%if (d2min <= R2max)
		%else
		%	
		%end
	end % interpolate_
	
	end % methods
end % classdef

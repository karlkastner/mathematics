% Sun Jul  6 11:24:55 WIB 2014
% Karl Kastner, Berlin
%
%% interpolator super-class
classdef Interpolator < handle
	properties
		% cut-off radius for source points around target points
		% this is transformed into an ellipse by aspect_ratio
		Rmax
		% cut-off radius around the verification points
		% Rmin < Rmax, e.g the cinclusion zone is a ring around the
		% verification points
		Rmin
		% order of the interpolation
		order
		% number of verification points
		nverify
%		ecov
		% tolerance for iterative solvers
		rtol
		% aspect_ratio converted to axis scale factors
		% subject to the condition (s(1)/s(2) = a and s(1)^2+s(2)^2 = 1)
		s
		qtree;
	end % properties
	methods
	% constructor
	function obj = Interpolator(Rmax,Rmin,order,nverify,aspect_ratio)
		obj.Rmax    = Rmax;
		obj.Rmin    = Rmin;
		obj.order   = order;
		obj.nverify = nverify;
		obj.rtol    = 1e-4;
		% aspect ratio anisotropic axes
		% r^2 = [dx dy][s1 0; 0 s2]^2[dx dy]' = s1^2 dx^2 + s2^2 dy^2
		% s1/s2 = aspect ratio
		% in transformed coordiantes : dx^2 + dy^2 = r^2
		% in untransformed coordinates : s1^2dx^2 + s2^2 dy^2 = r^2
		%                              : s1^2dx^2 + s2^2 (1-dx^2) = 1
		%                              : s1^2dx^2 + s1^2/a^2 (1-dx^2) = 1
		%			       : s1^2 = a^2/(a^2*dx^2 + (1-dx^2))
		% choose dx = 1/sqrt(2)        : s1^2 = 2a^2/(a^2 + 1)
		a = aspect_ratio;
		obj.s = sqrt([2*a*a/(1+a*a) 2/(1+a*a)]);
	end % constructor

	% initiates interpolation for each target point
	% actual interpolation is defined in sub-classes
	function [Vt Et Ot jd obj] = interpolate(obj,Xs,Vs,Xt)
		disp('Starting interpolation');
		obj.qtree = Qtree_(Xs(:,1),Xs(:,2));
		Vt = zeros(size(Xt,1),1);
		Et = zeros(size(Xt,1),1);
		Ot = zeros(size(Xt,1),1);	
		jd = cell(size(Xt,1),1);
		% for each target point
		timer = tic();
		t_old = 0;
		for idx=1:size(Xt,1)
			t = toc(timer);
			if (t-t_old > 10)
				fprintf(1,'Progress: %f%% %fs\n',100*idx/size(Xt,1),t);
				t_old = t;
			end % if
			x0 = Xt(idx,:);
			%func = @(dx,dy) cfunc(dx,dy,x0(1),x0(2),obj.ecov.param);
			%[Vt(idx) Et(idx) Ot(idx) jd{idx}] = obj.interpolate_(Xs,Vs,x0,0,func);
			% rmin is zero for interpolation
			[Vt(idx) Et(idx) Ot(idx) jd{idx}] = obj.interpolate_(Xs,Vs,x0,0);
		end % for each target point
	end % interpolate

	% initiates interpolation at verification points
	function [Ev vdx] = verify(obj,Xv,Vv,Xt,Xs,Vs)
		disp('Starting verification');
		% allocate memory
		Ev = zeros(obj.nverify,1);
		% for a random subset of points
		vdx = randi(size(Xv,1),obj.nverify,1);
		% determine the verification radius
		% TODO here actually the Mahalanobis distance should be used
	%	disp('Building nearest neighbour tree');
	%	qtree = Qtree_(Xs(:,1),Xs(:,2));
	%	[idnn dnn] = qtree.nearest_neighbour(Xt(:,1),Xt(:,1));
	%	TODO, this should be done in interpolate
		% find nearest neighbour in source for each point in target set
		[id dnn] = knnsearch(Xs,Xt,'K',1); 
		obj.Rmin = nanmedian(dnn);
		fprintf(1,'Verification radius is %f\n',obj.Rmin);
		% for each verification point
		%R2min = obj.Rmin*obj.Rmin;
		%R2max = obj.Rmax*obj.Rmax;
		t = tic();
		for idx=1:obj.nverify
			if (toc(t) > 10)
				fprintf('Progress %f\n',idx/obj.nverify);
				t = tic();
			end
			x0 = Xv(vdx(idx),:);
			func = @(dx,dy) cfunc(dx,dy,x0(1),x0(2),obj.ecov.param);
			vs_ver = obj.interpolate_(Xs,Vs,x0,obj.Rmin);
%			vs_ver = obj.interpolate_(Xs,Vs,x0,obj.Rmin,func);
%	function [Vt Et Ot jd obj] = interpolate(obj,Xs,Vs,Xt)
			% error is difference of interpolated value and measured value
			Ev(idx) = vs_ver - Vv(vdx(idx));
		end % for idx
	end % function verify

	% Mahalanobis distance
	function [D obj] = dist(obj,X1,X2)
		D(:,1) = X1(:,1) - X2(:,1);
		D(:,2) = X1(:,2) - X2(:,2);
		D(:,1) = obj.s(1)*D(:,1);
		D(:,2) = obj.s(2)*D(:,2);
	end

	end % methods
end % classdef Interpolator


% Sun 23 Jul 12:22:37 CEST 2017
%% class for regularized interpolation (Tikhonov) on a triangulation
%% (unstructured mesh)
classdef RegularizedInterpolator3 < handle
	properties
		% do not initialize the object here, matlab bug
		mesh;

		% values at mesh points
		vali = struct();
		
		% standard errors at mesh points
		msei = struct();

		% Thikonov regularisation factor
		lambda = [0 0 0];

		% order of Thikonov regularisation
		order = 1
		
		n = [];
		mseflag = false;
	end % properties
	methods
		% constructor
		function obj = RegularizedInterpolator3(varargin)
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
			obj.mesh = MMesh();
		end

		function obj = remesh(obj,xlim,ylim,zlim,n)
			% set up the mesh
			L   = [xlim(2)-xlim(1), ylim(2)-ylim(1), zlim(2)-zlim(1)];
			xy0 = [xlim(1),ylim(1),zlim(1)];
			obj.mesh = MMesh.generate_uniform_tetra(n,L,xy0);
			obj.n = n;
		end % remesh

		function [mse msei f2 obj] = interpolation_error(obj,name,varargin)
			if (nargin()<3)
				name = 'default';
			end
			[mse msei f2] = obj.mesh.interpolation_error_3d(obj.vali.(name)(:,1),varargin{:});
		end % interpolation_error

		function [I obj] = integrate(obj,name)
			if (nargin() < 2)
				name = 'default';
			end
			I = obj.mesh.integrate_3d(obj.vali.(name));
		end % integrate
	end % methods
end % RegularizedInterpolator3


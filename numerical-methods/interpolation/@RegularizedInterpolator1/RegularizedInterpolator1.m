% Thu 22 Jun 08:25:49 CEST 2017
%% class for regularized interpolation (Thikonov) on a 1D mesh
classdef RegularizedInterpolator1 < handle
	properties
		% This is a severe matlab bug, declaring mesh here
		% synchronises the property between all objects
		% mesh = MMesh();
		mesh

		% values at mesh points
		vali = struct();
		
		% standard errors at mesh points
		msei = struct();

		% Thikonov regularisation factor
		lambda = 0

		% order of Thikonov regularisation
		order = 1
		
		n = [];
	end % properties
	methods
		% constructor
		function obj = RegularizedInterpolator1(varargin)
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
			obj.mesh = MMesh();
		end

		function obj = remesh(obj,xlim,n)
			% set up the mesh
			%obj.mesh.generate_uniform_1d(n,L,xl);
			obj.mesh.generate_uniform_1d(n,xlim(2)-xlim(1),xlim(1));		
			obj.n = n;
		end

		function [mse, msei, f2, obj] = interpolation_error(obj,name,varargin)
			if (nargin()<2 || isempty(name))
				name = 'default';
			end
			[mse, msei, f2] = obj.mesh.interpolation_error_1d(obj.vali.(name)(:,1),varargin{:});
		end % interpolation_error

		function [I, obj] = integrate(obj,name)
			if (nargin() < 2)
				name = 'default';
			end
			I = obj.mesh.integrate_1d(obj.vali.(name));
		end % integrate

		function dx1 = dx1(obj)
			dx1 = range(obj.mesh.point(:,1))./(obj.n-1);
		end
	end % methods
end % RegularizedInterpolator1


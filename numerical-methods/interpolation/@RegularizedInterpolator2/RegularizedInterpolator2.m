% Sun 23 Jul 12:22:37 CEST 2017
%% class for regularized interpolation on an unstructures mesh (interpolation)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef RegularizedInterpolator2 < handle
	properties
		% do not initialize the object here, matlab bug
		mesh;

		% values at mesh points
		vali = struct();
		
		% standard errors at mesh points
		msei = struct();

		% Thikonov regularisation factor
		lambda = [0 0];

		% order of Thikonov regularisation
		order = 1
		
		n = [];
		mseflag = false;
	end % properties
	methods
		% constructor
		function obj = RegularizedInterpolator2(varargin)
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
			obj.mesh = UnstructuredMesh();
		end

		function obj = remesh(obj,xlim,ylim,n)
			% set up the mesh
			L = [xlim(2)-xlim(1), ylim(2)-ylim(1)];
			xy0 = [xlim(1),ylim(1)];
			obj.mesh = UnstructuredMesh.generate_uniform_triangulation(n,L,xy0);
			obj.n = n;
		end % remesh

		function [mse msei f2 obj] = interpolation_error(obj,name,varargin)
			if (nargin()<3)
				name = 'default';
			end
			[mse msei f2] = obj.mesh.interpolation_error_2d(obj.vali.(name)(:,1),varargin{:});
		end % interpolation_error

		function [I obj] = integrate(obj,name)
			if (nargin() < 2)
				name = 'default';
			end
			I = obj.mesh.integrate_2d(obj.vali.(name));
		end % integrate
	end % methods
end % RegularizedInterpolator2


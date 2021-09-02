% Wed 26 Aug 13:47:23 +08 2020
%
%% solve coupled first- and second-order 1D boundary-value problems
%
%
% discrete solver for boundary value problems
% employs the method of characteristics
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
%
classdef BVPS_Characteristic < handle
	properties
		% function returning ode-coefficients
		odefun
		% function returning boundary-coefficients
		bcfun
		% function returning initial-values
		inifun
	
		% struct containing solution and intermediate values
		% xc, dx
		out = struct()

		% number of channels
		nc
		% number of equations (frequency components)
		neq
		% 1d-coordinate of channel end points
		xi
		% number of grid points per channel
		nx
		% number of segments (nxc-1) per channel
		nxc
		% cumulated number of grid points
		ni
		% cumulated number of segments
		nci
		% cumulated number of segments for decomposed solution
		npi
		% cumulated number of segments for decomposed solution
		npii
		% order of equations
		oo

		% options
		opt  = struct();
		warn = struct('k',0,'kmax',3);

		% exponential function, can be replaced by @(x) (1+x)
		exp = @exp;

		% discretisation matrix buffer
		Abuf
		% associated counter
		nbuf
		% right hand side
		b

		% Ic
		% Dc

		% junction condition
		jfun
	end % properties
	methods
		function obj = BVPS_Characteristic()
		end % constructor
	end % methods	
end % classdef BVPS_Characteristic


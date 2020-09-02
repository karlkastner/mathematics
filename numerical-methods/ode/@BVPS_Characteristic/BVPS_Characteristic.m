% Wed 26 Aug 13:47:23 +08 2020
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
		odefun
		bcfun
		inifun
	
		% xc, dx,
		out = struct()

		nc
		xi
		nx
		nxc
		neq
		ni
		nci
		npi
		npii
		oo

		opt = struct();
		warn = struct('k',0,'kmax',3);

		exp = @exp;

		% discretisation matrix buffer
		Abuf
		nbuf
		b

		jfun
	end % properties
	methods
		function obj = BVPS_Characteristic()
		end % constructor
	end % methods	
end % classdef BVPS_Characteristic


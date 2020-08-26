% Wed 26 Aug 13:47:23 +08 2020
%
% discrete solver for boundary value problems
% employs the method of characteristics
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

		exp = @exp;

		jfun
	end % properties
	methods
		function obj = BVPS_Characteristic()
		end % constructor
	end % methods	
end % classdef BVPS_Characteristic


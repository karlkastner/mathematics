% Mon 10 Jan 11:11:55 CET 2022
%
% class for analysis of remotely sensed and model generated vegetation patterns
%
classdef Spatial_Pattern < handle
	properties
		% pattern in real space
		b

		% mask of pattern in real space
		msk
		%b_hp;
		
		% domain size in m	
		L

		% spectral density
		S = struct();

		% autorrelatin
		R = struct();
		
		% real axes (m)
		x
		y
		r

		% frequency axes 1/m
		f = struct('x',[],'y',[],'r',[]);

		% windows
		w 

		% analysis options
		opt = struct( 'correctlocalmean', false ...
			     ,'fitmethod', 'ls'	...
			     ,'nb', [] ...
			     ... % levels for confidence intervals for density plots
			     , 'pp', [0.16,0.84] ...
			     , 'normalize', true ...
			     , 'r2_min', 0.5 ...
			     , 'nf_radial', 3 ...
			     ... % confidence level for periodicity
			     , 'confidence_level', 0.05 ...
			     , 'template', 'bartlett' ...
			     , 'fminscale', 1/8 ...
			     , 'fmaxscale', 4 ...
			     , 'mseg_scale', 20 ... % 6 % TODO, this can become critical, when pattern too short
		);
		% analysis results
		stat = struct( );
	end % properties
	methods
		function obj = Spatial_Pattern(varargin)
			for idx=1:2:length(varargin)
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end % for idx
		end % constructor Spatial_Pattern
		function n = n(obj)
			if (isvector(obj.b))
				n = length(obj.b);
			else
				n = size(obj.b);
			end
		end
		function df  = df(obj)
			% frequency interval (spacing of frequency bins)
			df = 1./obj.L;
		end
	end % methods
end % classdef Spatial_Pattern

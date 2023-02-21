% Mon 10 Jan 11:11:55 CET 2022
% Karl Kastner, Berlin
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% class for analysis of remotely sensed and model generated vegetation patterns
%
classdef Spatial_Pattern < handle
	properties
		% pattern in real space
		b
		b_square

		% mask of pattern in real space
		msk = struct('b',[],'f',[],'rot',[]);
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
			     , 'ns', 100 ...
			     , 'n_mem_blk', 2.5e8 ...
			     , 'objective', 'mise-cramer' ...
			     , 'xlim', [0, 2.5] ...
			     , 'pfmsk', 0.8 ...
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
		function lambda_c = lambda_c(obj)
			if (obj.stat.isisotropic)
				Sc = obj.stat.Sc.radial.clip;
				lambda_c = 1./obj.stat.fc.radial.clip;	
			else
				Sc = obj.stat.Sc.x.clip;
				lambda_c = 1./obj.stat.fc.x.clip;
			end % else of isiso
		end % lambda_c
	end % methods
end % classdef Spatial_Pattern

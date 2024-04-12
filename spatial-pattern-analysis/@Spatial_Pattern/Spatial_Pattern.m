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
		% pattern in the spatial domain
		b
		% extended to square domain by padding the mean value
		b_square

		% mask of pattern in real space
		msk = struct('b',[],'f',[],'rot',[]);
		
		% spatial extend (domain size) in m	
		L

		% periodogram and spectral density
		S = struct();

		% cumulative spectral density
		C = struct();

		% correlogram and autocorrelation
		R = struct();
		
		% axes int the spatial domain (m)
		x
		y
		r

		% axes in the frequency domain 1/m
		f = struct('x',[],'y',[],'r',[]);

		% windows for suppressing spurious frequency components
		w 

		% analysis options
		opt = struct( ... 
			       'datatype', 'single' ...
			     , 'n_mc', 100 ... % number of repetitions for monte-carlo
			     , 'nf_test_min', 3 ...
			     , 'nf_test_scale', 1 ...
			     , 'objective', 'least-squares' ...
			     , 'p_fmsk', 0.8 ...
			     , 'significance_level_a1', 0.05 ...
			     , 'test_for_periodicity', true ...
			     , 'xlim', [0, 2.5] ...
			     , 'transect', struct('mseg_scale', 20 ... % 6 % TODO, this can become critical, when pattern too short
						  , 'fitmethod', 'ls' ...
			                          ,'nb', [] ...
			     			  ... % levels for confidence intervals for density plots
					          , 'pp', [0.16,0.84] ...
			                          , 'normalize', true ...
			                          , 'r2_min', 0.5 ...
			                          , 'template', 'bartlett' ...
			                          , 'fminscale', 1/8 ...
			                          , 'fmaxscale', 4 ...
				           ) ...
			    , 'scalefield', 'hp' ...
			    , 'weight', true ...
		);

		% analysis statistics
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
		function Sc = Sc(obj)
			if (obj.stat.isisotropic)
				Sc = obj.stat.Sc.radial.(obj.opt.scalefield);
			else
				Sc = obj.stat.Sc.x.(obj.opt.scalefield);
			end % else of isiso
		end
		function lambda_c = lambda_c(obj)
			if (obj.stat.isisotropic)
				lambda_c = 1./obj.stat.fc.radial.(obj.opt.scalefield);	
			else
				lambda_c = 1./obj.stat.fc.x.(obj.opt.scalefield);
			end % else of isiso
		end % lambda_c
	end % methods
end % classdef Spatial_Pattern

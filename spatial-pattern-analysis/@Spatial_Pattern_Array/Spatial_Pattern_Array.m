% 2022-09-26 14:32:10.449630621 +0200
% Karl Kästner, Berlin
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
%
%% container class for Spatial_Pattern objects
classdef Spatial_Pattern_Array < handle
	properties

		% array of spatial pattern objets
		sp_a	   = Spatial_Pattern();
		% centroid coordinates
		centroid   = [];
		% sampling resolution
		dx_sample
		fr_max_target;
		% processing errors
		error_C    = {};
		% file names
		filename_C = [];
		region_id  = [];
		region_C   = '';
		runtime    = [];
		type       = '';
		%e(end:n)           = {[]};
		%filename_C(end:n)   = {''};
		opt = struct( ... 
			      'test_for_periodicity',true ...
			     ,'area_max', 4e6 ...
			     ,'cmd_str', 'LD_LIBRARY_PATH= python3 vegetation_patterns_fetch_polygons_google_satellite.py %s %g' ...
			     ,'dx0', 2 ...
			     ,'dx_max', 4 ...
			     ,'dx_min', 0.5 ...
			     ,'regularity_min', 0.4 ... % minimum regularity of patterns included in the analysis
			     ,'reload', true ...
			     ,'tail_trim_scale', 3 ...
			     ,'wavelength_max', 500 ...% maximum wavelength for a pattern to be included
			     ,'wavelength_min', 1 ...% minimum wavelength for a pattern to be included
			     , 'skip', 1 ...
			    );
		imgbase  = 'google-satellite_';
	end

	methods
	function obj = Spatial_Pattern_Array(centroid,dx_sample)
		if (nargin()>0)
			obj.centroid = centroid;
		end
		if (nargin()>1)
			obj.dx_sample = dx_sample;
		end
		obj.runtime = NaN(obj.n,1);
		obj.region_id = ones(obj.n,1);
		obj.error_C = cell(obj.n,1);
		obj.filename_C = repmat({''},obj.n,1);

	end

	% pseudo members
	function n = n(obj)
		n = size(obj.centroid,1);
	end
	function y = area_msk(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.area_msk',NaN)),obj.sp_a(idx));
	end
	function y = contrast(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.contrast',NaN)),obj.sp_a(idx));
	end
	function y = coverage(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.coverage',NaN)),obj.sp_a(idx));
	end
	function y = date(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.date',NaN)),obj.sp_a(idx));
	end
	function y = isisotropic(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(double(getfield_try(x,'stat.isisotropic',NaN))),obj.sp_a(idx));
	end
	function y = intS_hp_sig(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.stati.intS_hp_sig',NaN)),obj.sp_a(idx));
	end
	function y = L_eff_r(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.L_eff.r',NaN)),obj.sp_a(idx));
	end
	function y = L_eff_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.L_eff.x',NaN)),obj.sp_a(idx));
	end
	function y = L_eff_y(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.L_eff.y',NaN)),obj.sp_a(idx));
	end

	function y = p_periodic(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.p_periodic',NaN)),obj.sp_a(idx));
	end
	function y = p_S_hp(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.p_S_hp',NaN)),obj.sp_a(idx));
	end
	function y = wavelength_r(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) 1./double(getfield_try(x,'stat.fc.radial.hp',NaN)),obj.sp_a(idx));
	end
	function y = wavelength_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) 1./double(getfield_try(x,'stat.fc.x.hp',NaN)),obj.sp_a(idx));
	end
	function y = Sc_r(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.Sc.radial.hp',NaN)),obj.sp_a(idx));
	end
	function y = Sc_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.Sc.x.hp',NaN)),obj.sp_a(idx));
	end
	function y = Sc_y(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,'stat.Sc.y.hp',NaN)),obj.sp_a(idx));
	end
	function y = regularity_r(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Sc_r(idx)./obj.wavelength_r(idx);
	end
	function y = regularity_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Sc_x(idx)./obj.wavelength_x(idx);
	end
	function y = regularity_y(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Sc_y(idx)./obj.wavelength_y(idx);
	end
	end % methods
end % class Spatial_Pattern_Array


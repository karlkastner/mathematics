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
		%dx_sample
		level_a
		fr_max_target;
		% processing errors
		error_C    = {};
		% file names
		filename_C = [];
		% world region indices
		region_id  = [];
		% world region strings
		region_C   = '';
		% runtime of analyzis
		runtime    = [];
		% either isotropic or anisotropic
		type       = '';
		% object for knnsearch for matching input files based on centroid coordinates
		knnobj_C   = [];
		% pathes and filenames of satellite images
		img_C      = [];
		% options
		opt = struct( ...
			       'analyze', true ... % (new) satellite images are only anlyzed if true
			     , 'area_max', 4e6 ... % maximum area in pixels, to avoid out of memory
			     , 'base_str', ''  ... % base of input/output files 
			     , 'cmd_str', 'LD_LIBRARY_PATH= python3 vegetation_patterns_fetch_polygons_google_satellite.py %s %g' ...
			     , 'd_max', 1e-4 ...   % maximum distance for matching files based on centroid coordinates
			     , 'dx0', 2 ...	   % initial sampling resolution
			     , 'dx_max', 4 ...	   % max sampling resolution
			     , 'dx_min', 0.5 ...   % min sampling resolution
			     , 'field', 'con' ...  % state analysis results based on consistently estimated densities
			     , 'folder_img', 'img/%s/%g/' ... % folder for satellite images
			     , 'folder_mat', 'img/%s/analysis/%g/' ... % folder for individual analysis results of each image
			     , 'imgbase', 'google-satellite_' ... % basename for satellite images
			     , 'match_by_knn', true ... % if true, file names are matched by nn search of centroid coordinates
			     , 'regularity_min', 0.4 ... % minimum regularity of patterns included in the analysis
			     , 'reload', true ... % reuse anlysis computed during previous run
			     , 'skip', 1 ...      % skip every skip-1 pattern for quick test
			     , 'tail_trim_scale', 3 ... % skale for heuristically determining the sampling resolution
			     , 'test_for_periodicity', true ... % test for periodicity during analysis
			     , 'wavelength_max', 500 ... % maximum wavelength for a pattern to be included
			     , 'wavelength_min', 1 ... % minimum wavelength for a pattern to be included
			    );
	end % properties

	methods
	function obj = Spatial_Pattern_Array(centroid,level_a,dx_sample)
		%dx_sample)
		if (nargin()>0)
			obj.centroid = centroid;
		end
		if (nargin()>1 && ~isempty(level_a))
			obj.level_a = level_a;
		end
		if (nargin()>2)
			obj.level_a = cvec(obj.resolution2level(dx_sample));
		end
		obj.runtime = NaN(obj.n,1);
		obj.region_id = ones(obj.n,1);
		obj.error_C = cell(obj.n,1);
		obj.filename_C = repmat({''},obj.n,1);

	end
	% pseudo members
	function level = resolution2level(obj,dx)
		level = round(-log2(dx/obj.opt.dx0))+1;
		level(dx <= 0) = NaN;
	end
	function level_min = level_min(obj)
		level_min = obj.resolution2level(obj.opt.dx_max);
		%level_min = round(-log2(obj.opt.dx_min/obj.opt.dx0))+1;
	end
	function level_max = level_max(obj)
		level_max = obj.resolution2level(obj.opt.dx_min);
		%level_max = round(-log2(obj.opt.dx_max/obj.opt.dx0))+1;
	end
	function dx = resolution(obj,level)
		dx = obj.opt.dx0*2.^(-level+1);
	end
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
		y = arrayfun(@(x) getfield_try(x,'stat.date',''),obj.sp_a(idx),'uniformoutput',false);
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
		y = arrayfun(@(x) 1./double(getfield_try(x,['stat.fc.radial.',obj.opt.field],NaN)),obj.sp_a(idx));
	end
	function y = wavelength_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) 1./double(getfield_try(x,['stat.fc.x.',obj.opt.field],NaN)),obj.sp_a(idx));
	end
	function y = Src(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,['stat.Sc.radial.',obj.opt.field],NaN)),obj.sp_a(idx));
	end
	function y = Sxpc(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,['stat.Sc.xp.',obj.opt.field],NaN)),obj.sp_a(idx));
	end
	function y = Syc(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = arrayfun(@(x) double(getfield_try(x,['stat.Sc.y.',obj.opt.field],NaN)),obj.sp_a(idx));
	end
	function y = regularity_r(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Src(idx)./obj.wavelength_r(idx);
	end
	function y = regularity_x(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Sxpc(idx)./obj.wavelength_x(idx);
	end
	function y = regularity_y(obj,idx)
		if (nargin()<2) idx=1:obj.n; end
		y = obj.Syc(idx)./obj.wavelength_y(idx);
	end
	end % methods
end % class Spatial_Pattern_Array


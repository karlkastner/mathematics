% Wed 18 May 13:49:09 CEST 2022
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
%% read an image of a pattern, its mask and geospatial data
function [g, alpha, obj] = imread(obj,filename)
	g = GeoImg();
	g.read(filename);

	b = g.img;

	%[img, map, alpha] = imread(filename);
	if (~isempty(g.map))
		error('Spatial_Pattern:ImageIsIndexed','Convert indexed images to grayscale or RGB first');
	end
	switch (obj.opt.datatype)
	case {'single'}
		b = single(b);
	case {'double'}
		b = double(b);
	otherwise
		error('Spatial_Pattern:DataType','Dataype must be single or double');
	end
	if (3 == ndims(b))
		% sloppy conversion to grayscale
		b = mean(b,3);
	end

	% invert to get biomass proxy
	b = max(b,[],'all')-b;
	% stretch contrast to 1
	b = b./max(b,[],'all');

	n = size(b);

	obj.stat.pgw = g.pgw;
	obj.stat.xy0 = g.xy0;
	obj.stat.dxy = g.dxy;
	obj.stat.bangle_rad = g.angle;
	obj.stat.date = g.info.FileModDate;

	try
		% try to load mask from file
		mskfilename = [filename(1:end-4),'_mask.tif'];
		msk = imread(mskfilename);
		nmsk = size(msk);
		if (rms(nmsk(1:2)-n(1:2))>0)
			mskfilename
			error('SpatialPattern:MaskSize','Mask size must match image size');
		end
	catch
		% use alpha data
		if (~isempty(g.alpha))
			% note that we can even alias the maks, but we do not do this here
			msk = g.alpha > 0;
		else
		msk = true(n);
		end
	end % catch of try

	obj.L   = n.*obj.stat.dxy;
	obj.b   = b;
	obj.msk.b = msk;
	obj.stat.filename = filename;
end % imread


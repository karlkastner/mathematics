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
%% read an image file containing a pattern, mask and geospatial data
function [img, alpha, obj] = imread(obj,filename)
	g = GeoImg();
	g.read(filename);
	%[img, map, alpha] = imread(filename);
	if (~isempty(g.map))
		error('Convert indexed images to grayscale or RGB first');
	end
	switch (obj.opt.datatype)
	case {'single'}
		img = single(g.img);
	case {'double'}
		img = double(g.img);
	otherwise
		error('Dataype must be single or double');
	end
	if (3 == ndims(img))
		% sloppy conversion to grayscale
		img = mean(img,3);
	end

	% invert to get biomass proxy
	b = max(img,[],'all')-img;
	% stretch contrast to 1
	b = b./max(b,[],'all');

	n = size(b);

	%try
		%out = readpgw([filename(1:end-4),'.pgw']);
		%dxy = g.dxy;
		obj.stat.pgw = g.pgw;
		obj.stat.xy0 = g.xy0;
		obj.stat.dxy = g.dxy;
		obj.stat.bangle_rad = g.angle;
		obj.stat.date = g.info.FileModDate;
	%catch
	%	dxy = 1;
	%	obj.stat.pgw = [];
	%	obj.stat.dxy = NaN;
	%	obj.stat.xy0 = [NaN,NaN];
	%	obj.stat.bangle_rad = NaN;
	%end

	try
		% try to load mask from file
		mskfilename = [filename(1:end-4),'_mask.tif'];
		msk = imread(mskfilename);
		nmsk = size(msk);
		if (rms(nmsk(1:2)-n(1:2))>0)
			mskfilename
			disp('mask size does not match');
			error('msksize');
		end
	catch
		% use alpha data
		%disp('Unable to load mask, analyzing whole area');
		if (~isempty(g.alpha))
			% note that we can even alias the maks, but we do not do this here
			msk = g.alpha > 0;
		else
		msk = true(n);
		end
	end % catch of try
	%end % else of ~ isempty alpha

	obj.L   = n.*obj.stat.dxy;
	obj.b   = b;
	obj.msk.b = msk;
end


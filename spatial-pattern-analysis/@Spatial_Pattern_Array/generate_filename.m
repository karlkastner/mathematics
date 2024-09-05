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
% TODO store d in output
function [folder_img,filename,folder_a] = generate_filename(obj,idx,level)
	resolution = obj.resolution(level);
	lid = level+1;

	folder_img = sprintf(obj.opt.folder_img,obj.type,resolution); 
	folder_a   = sprintf(obj.opt.folder_mat,obj.type,resolution); 
	%folder     = ['img/',obj.type,'/',num2str(resolution_),'/'];
	%folder_a   = ['img/',obj.type,'/analysis/',num2str(resolution_),'/'];
	x_str      = sprintf('%+09.5f',obj.centroid(idx,2));
	y_str      = sprintf('%+010.5f',obj.centroid(idx,1));
	if (obj.opt.match_by_knn)
	% match pattern image file by coordinates
	% facilitated by nearest neibour search
	% the pattern image files name specifies the coordinates, unfortunately 
	% the coordinates computed in python and matlab differ slightly due to
	% round-off error, due to rounding of 0.000 or 9.999 this can affect all
	% digits in the worst case, therefore files are matched here by a nearest
	% neighbour search
	if (length(obj.knnobj_C) < lid || isempty(obj.knnobj_C{lid}))
		printf('Creating search object for the level %g, dx = %g\n',level,resolution);
		obj.img_C{lid} = dir([folder_img,'/google-satellite_*png']);
		C = [arrayfun(@(x) str2num(x.name(18:26)), obj.img_C{lid}), ...
	             arrayfun(@(x) str2num(x.name(28:37)), obj.img_C{lid})];
		obj.knnobj_C{lid} = createns(C);
	end
	[id,d] = knnsearch(obj.knnobj_C{lid},obj.centroid(idx,2:-1:1));
	if (d < obj.opt.d_max)
		filepath = [obj.img_C{lid}(id).folder,filesep,obj.img_C{lid}(id).name];
	else
		% return default filename	
		xy_str = [x_str,'_',y_str];
		filename = [obj.opt.imgbase,xy_str,'.png'];
		filepath = [folder_img,filesep(),filename];
	end
	else % match filenames directly
	try
		% try wildcards
		xy_str = [x_str(1:end-1),'*_',y_str(1:end-1),'*'];
		filename = [obj.opt.imgbase,xy_str,'.png'];
		filepath = chomp(ls([folder_img,filesep(),filename]));
	catch e
		% return default filename	
		xy_str = [x_str,'_',y_str];
		filename = [obj.opt.imgbase,xy_str,'.png'];
		filepath = [folder_img,filesep(),filename];
	end
	end % match dy direct name
	filename = basename(filepath);
end


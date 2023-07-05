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
% TODO nearest neighbour matching
function [folder,filename,folder_a] = generate_filename(obj,idx,resolution_)
	folder     = ['img/',obj.type,'/',num2str(resolution_),'/'];
	folder_a   = ['img/',obj.type,'/analysis/',num2str(resolution_),'/'];
	x_str      = sprintf('%+09.5f',obj.centroid(idx,2));
	y_str      = sprintf('%+010.5f',obj.centroid(idx,1));
	try
		% try wildcards
		xy_str = [x_str(1:end-1),'*_',y_str(1:end-1),'*'];
		filename = [obj.imgbase,xy_str,'.png'];
		filepath = chomp(ls([folder,filesep,filename]));
	catch e
		% return default filename	
		xy_str = [x_str,'_',y_str];
		filename = [obj.imgbase,xy_str,'.png'];
		filepath = [folder,filesep,filename];
	end
	filename = basename(filepath);
end


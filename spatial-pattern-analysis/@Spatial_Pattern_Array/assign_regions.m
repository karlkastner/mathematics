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
function assign_regions(obj,filename)
	regions_shp = Shp.read(filename);

	% assign regions
	region_C = {regions_shp.region};
	obj.region_C   = vertcat('other',region_C(:)); 
	% other has index 1
	obj.region_id  = 1 + Shp.inpolygon(regions_shp,obj.centroid(:,1),obj.centroid(:,2));
end
	

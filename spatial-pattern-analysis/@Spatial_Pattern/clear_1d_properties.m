% 2022-09-26 14:32:10.449630621 +0200
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
function obj = clear_1d_properties(obj)
		% clear 2d properties first
		obj.clear_2d_properties();
		% remove 2d grid fields to save memory and disk space
		obj.x = [];
		obj.y = [];
		obj.r = [];
		obj.S = [];
		obj.C = [];
		obj.R = [];
		obj.f = []; % x,y,r,a
end


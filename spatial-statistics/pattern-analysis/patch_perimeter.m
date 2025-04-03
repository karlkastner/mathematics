% 2024-11-28 14:00:25.890335569 +0100
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
function pp = patch_perimeter(z,val)
	pp = (z==val) ...
		& ( left(z,1) ~= val ...
                       | right(z,1) ~= val ...
                       | up(z,1) ~= val ...
                       | down(z,1) ~= val ...
		   );
end

% 2025-03-11 14:22:36.530277810 +0100
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
%
function S = resample_density_by_averaging_2d(S,no)
	if (length(no)<2)
		no(2) = no(1);
	end
	% for first and second dimension
	for idx=1:2
		S = resample_density_by_averaging(S,no(idx),idx);
	end
end


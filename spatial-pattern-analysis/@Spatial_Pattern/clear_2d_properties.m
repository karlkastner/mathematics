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
function obj = clear_2d_properties(obj)
	% remove 2d grid fields to save memory and disk space
	obj.b = [];
	obj.b_square = [];
	obj.msk      = [];
	obj.stat.qq = [];
	for field = {'hat','hp','bar','con'}
		obj.S.(field{1}) = [];
		obj.R.(field{1}) = [];
		obj.S.rot.(field{1}) = [];
		obj.R.rot.(field{1}) = [];
	end % for field
	obj.f.rr = [];
	obj.f.tt = [];
	%obj.w    = [];
	obj.stat.stati.p1_all = [];
	obj.stat.stati.pn_all = [];
end


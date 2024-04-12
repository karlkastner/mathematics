% 2024-01-22 10:01:19.374305026 +0100
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
function [x,y,rr] = real_axis_2d(n,L);
	x = real_axis_1d(n(1),L(1));
	y = real_axis_1d(n(2),L(2));
	if (nargout()>2)
		rr = hypot(x,y');
	end
end

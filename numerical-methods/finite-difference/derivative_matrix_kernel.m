% Wed 22 Nov 12:13:21 CET 2023
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
%
function d = derivativew_matrix_kernel(l,d,x0)
	x  = -(l-1)/2:(l-1)/2;
	v  = vander_1d(x,l-1);
	vd = vanderd_1d(x0,l-1,d);
	d  = (vd*inv(v));
end


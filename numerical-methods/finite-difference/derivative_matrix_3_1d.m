% Wed 22 Nov 09:36:22 CET 2023
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
function D = derivative_matrix_3_1d(L,n)
	dx = L/(n-1);
	s = 1/dx^3;
	d = s*[-0.5    1.0         0   -1.0    0.5];
	D = spdiags(ones(n,1)*d,-2:2,n,n);
	% extrapolation at boundary
	D(1,1:5)           = s*[ -2.5    9.0  -12.0    7.0   -1.5];
	D(2,1:5)           = s*[ -1.5    5.0   -6.0    3.0   -0.5];
	D(end-1,end-4:end) = s*[  0.5   -3.0    6.0   -5.0    1.5];
	D(end,end-4:end)   = s*[  1.5   -7.0   12.0   -9.0    2.5];
end


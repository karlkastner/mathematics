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
function D = derivative_matrix_4_1d(L,n)
	dx = L/(n-1);
	s = 1/dx^4;
	d = s*[1.0   -4.0    6.0   -4.0    1.0];
	D = spdiags(ones(n,1)*d,-2:2,n,n);
	% extrapolation at boundary, constant
	% linear extrapolation would require a 6 point kernel
	D(1,1:5)           = s*d;
	D(2,1:5)           = s*d;
	D(end-1,end-4:end) = s*d;
	D(end,end-4:end)   = s*d;
end


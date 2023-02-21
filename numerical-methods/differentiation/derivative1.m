% Fri 15 Dec 11:00:30 CET 2017
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
%
%
%% first derivative on variable mesh
%% second order accurate
% function dy_dx = derivative1(x,y,order)
function dy_dx = derivative1(x,y,order)
	istransposed = false;
	if (isvector(y) && isrow(y))
		y = cvec(y);
		istransposed = true;
	end

	D1 = derivative_matrix_1_1d(x,[]);

	dy_dx = D1*y;

	if (istransposed)
		dy_dx = dy_dx.';
	end
end


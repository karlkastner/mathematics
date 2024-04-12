% 2024-01-04 20:22:25.391693409 +0100
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
function [R1,R2t,R1_,R2_] = downsampling_matrix_2d(n)
	R1  = downsampling_matrix_1d(n(1),'mg');
	R2t = downsampling_matrix_1d(n(2),'mg')';
	if (nargout()>2)
		I1 = speye(n(1))
		I2 = speye(n(2))
		R1_ = kron(R1,I2);
		R2_ = R1_';
		%R2_ = kron(I1,R2t);
	end
end


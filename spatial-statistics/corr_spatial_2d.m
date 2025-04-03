% Fri 31 Jan 12:04:54 CET 2025
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
% spatial correlation in two dimensions
% assumes circular bc
function rho = corr_spatial_2d(z)
	% zij = 1/4*(zl + zr + zu + zd)*rho + eps
	% ((1 - r/4)I - r/4*D) z = eps
	rho = zeros(size(z,3),1);
	for idx=1:size(z,3)
		rhs      = (   up(z(:,:,idx),1) + down(z(:,:,idx),1) ...
			     + left(z(:,:,idx),1) + right(z(:,:,idx),1) ...
			   );
		rho(idx) = (rhs(:) \ flat(z(:,:,idx)))/4;
	end
end



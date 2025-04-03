% Fri 31 Jan 11:52:20 CET 2025
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
function rho = corr_spatial_1d(z)
	% zij = 1/2*(zl + zr)*rho + eps
	% ((1 - r/2)I - r/2*D) z = eps
	rho = zeros(size(z,2),1);
	z = z-mean(z);
	for idx=1:size(z,2)
		rhs      = (up(z(:,idx),1) + down(z(:,idx),1));
		%rho(idx,1) = (z(:,idx) \ rhs(:))/2;
		rho(idx,2) = (z(:,idx)'*rhs)/(z(:,idx)'*z(:,idx))/2;
	end
end


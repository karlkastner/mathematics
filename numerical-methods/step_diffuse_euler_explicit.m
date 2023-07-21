% Tue 18 Jul 09:50:19 CEST 2023
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
function z = step_diffuse_euler_explicit(dt,z,dx,e)
% TODO 2D
	n  = length(z);
	D2 = derivative_matrix_2_1d(n,dx,2,'circular','circular',true);
	z  = z + dt*e*D2*z; 
end



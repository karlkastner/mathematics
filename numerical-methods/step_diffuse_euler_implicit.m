% Tue 18 Jul 09:49:29 CEST 2023
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
% function z = diffuse_euler_implicit(dt,z,n,dx,e)
function z = diffuse_euler_implicit(dt,z,n,dx,e)
	Ix = speye(n(1));
	e_D2x = (dt*e(1))*derivative_matrix_2_1d(n(1),dx(1),2,'circular','circular',true);
	if (1 == length(n))
		A = (Ix - e_D2x);
	else
		Iy    = speye(n(2));
		e_D2y = (dt*e(2))*derivative_matrix_2_1d(n(2),dx(2),2,'circular','circular',true);

		I = speye(prod(n)) 
		e_D2 = kron(e_D2x,Iy) + kron(Ix,e_D2y);

		A = I - e_D2;
	end
	z = A \ z; 
end


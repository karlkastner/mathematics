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
function z = advect_euler_explicit(dt,z,n,dx,a)
	Ix = speye(n(1));
	D1x = (dt*a(1))*derivative_matrix_1_1d(n(1),dx(1),sign(a(1)),'circular','circular',true);
	if (1 == length(n))
		A = (Ix -a_D1x);
	else
		Iy  = speye(n(2));
		a_D1y = (dt*a(2))*derivative_matrix_2_1d(n(2),dx(2),sign(a(2)),'circular','circular',true);

		a_D1 = kron(I_y,a_D1x) + kron(a_D1y,Ix);

		A = speye(prod(n)) - a_D1;
		
		%A = kron(Ix - dt*e(1)*D2x, Iy) + kron(Ix,Iy - dt*e(2)*D2y);
	end
	z = A \ z; 
end


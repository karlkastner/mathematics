% Tue 18 Jul 09:52:07 CEST 2023
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
function z = diffuse_trapezoidal(dt,z,n,dx,e)
	% TODO use krylov
	Ix  = speye(n(1));
	D2x = derivative_matrix_2_1d(n(1),dx(1),2,'circular','circular',true);
	if (1 == length(n))
		A_l = Ix - 0.5*dt*e(1)*D2x;
		A_r = Ix + 0.5*dt*e(1)*D2x;
	else
		Iy  = speye(n(2));
		D2y = derivative_matrix_2_1d(n(2),dx(2),2,'circular','circular',true);

		A_l = kron(Ix - 0.5*dt*e(1)*D2x, Iy) + kron(Ix, Iy - 0.5*dt*e(2)*D2y);
		A_r = kron(Ix + 0.5*dt*e(1)*D2x, Iy) + kron(Ix, Iy + 0.5*dt*e(2)*D2y);
	end
	z = A_l \ (A_r*z);
end


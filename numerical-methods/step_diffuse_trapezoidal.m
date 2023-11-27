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
function [z,A_l,A_r] = step_diffuse_trapezoidal(dt,z,n,dx,e)
	I  = speye(n(1));
	D2 = (0.5*dt*e(1))*derivative_matrix_2_1d(n(1),dx(1),2,'circular','circular',true);
	if (2 == length(n))
		Iy  = speye(n(2));
		D2y = (0.5*dt*e(2))*derivative_matrix_2_1d(n(2),dx(2),2,'circular','circular',true);

		D2 = kron(D2y,I) + kron(Iy,D2);

		I = speye(prod(n));
	end
	A_l = (I - D2);
	A_r = (I + D2);
	z   = A_l \ (A_r*z);
end


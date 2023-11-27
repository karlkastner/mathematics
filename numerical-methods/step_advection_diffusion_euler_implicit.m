% Wed 19 Jul 17:38:03 CEST 2023
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
function [z,A] = step_advection_diffusion_implicit_euler(dt,dx,n,z,a,e)
	% TODO use krylov
	D1 = (dt*a(1))*derivative_matrix_1_1d(n(1),dx(1),2,'circular','circular',true);
	D2 = (dt*e(1))*derivative_matrix_2_1d(n(1),dx(1),2,'circular','circular',true);
	I  = speye(n(1));
	if (length(n)>1)
		D1y = (0.5*dt*a(2))*derivative_matrix_1_1d(n(2),dx(2),2,'circular','circular',true);
		D2y = (0.5*dt*e(2))*derivative_matrix_2_1d(n(2),dx(2),2,'circular','circular',true);
		Iy  = speye(n(2));
		D1  = kron(Iy,D1) + kron(D1y,I);
		D2  = kron(Iy,D2) + kron(D2y,I);
		I   = speye(prod(n));
	end
	A = (I - D1 - D2);
	z = A \ z;
end


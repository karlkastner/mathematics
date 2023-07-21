% 2023-07-17 10:16:00.698550380 +0200
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
% semi-analytic euler forward step for the inhomogeneous reaction equation
function z = step_euler_forward_exp(t,z,dt,afun,bfun)
	% semi-analytic
	a = afun(t,z);
	b = bfun(t,z);
	%y = exp(a*dt)*y0 + (exp(a*dx) - 1)*b/a
	%y = exp(a*dt)*y0 + (exp(a*dx) - 1)*b/a
	z = exp(dt*a).*z + dt*b;
end


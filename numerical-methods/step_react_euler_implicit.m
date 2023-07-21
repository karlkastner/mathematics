% 2023-07-17 10:15:25.081983737 +0200
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
function z = step_react_euler_implicit(t,z,dt,afun,bfun)
% implicit euler
		error('TODO interate if necessary');
		a = afun(t,z);
		b = bfun(t,z);
		z = (z + dt*b) ./ (1 - dt*a);
%		if (max(abs(a*dt))>1)
%			warning([num2str(max(abs(a*dt))),'max(|a dt|)']);
%		end
end

% 2023-07-18 20:19:19.600831047 +0200
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

n = 10*[1,1]; dx=[1,1]; dt=1; v=[0.5,0.5]; a = advection_kernel(dx,dt,n,v); [[NaN;(0:n-1)'],[0:n-1;a]]

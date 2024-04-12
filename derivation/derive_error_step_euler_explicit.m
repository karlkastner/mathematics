% Tue 14 Nov 11:42:17 CET 2023
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

% dy/dt = y' = -a y
% euler implicit
% y(t+dt) = y(t) + dt*y'(t+dt) 

% yi+1 - yi = dt(y'i+1)
% y(t+dt) = y(t) + y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y''' dt^3 
% y'(t+dt) = y'(t) + y''(t)*dt + 1/2 y''' dt^3
% y(t) + y'(t+1)*dt + E = y(t) + y'(t)*dt + 1/2*y''(t)*dt^2 
%    0 + (y'(t) + y''(t)*dt + 1/2 y''' dt^2)*dt + E =  0 + y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y''' dt^3
%    0 + ( 0 + y''(t)*dt + 1/2 y''' dt^2)*dt + E =  0 + 0 + 1/2*y''(t)*dt^2 + 1/6 y''' dt^3
%    E =  - 1/2*y''(t)*dt^2 - 1/3 y'''*dt^3

syms y d1y d2y d3y d4y dt d5y e

ynext  = y  + d1y*dt + 1/2*d2y*dt^2 + 1/6*d3y*dt^3 + 1/24*d4y*dt^4 + 1/120*d5y*dt^5
%dynext = d1y + d2y*dt + 1/2*d3y*dt^2 + 1/6*d4y*dt^3 + 1/24*d5y*dt^4
ynum = y + dt*d1y

eq = ynext == ynum + e
solve(eq,e)

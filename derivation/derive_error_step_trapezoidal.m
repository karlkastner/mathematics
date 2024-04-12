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

% dy/dt = y'
% trapezoidal
%  y(t+dt) = y(t) + dt/2*(y'(t)+y'(t+dt))

% yi+1 - yi = dt/2(y'i+1 + y'i)
% y(t+dt) = y(t) + y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y'''(t)*dt^3
% y'(t+dt) = y'(t) + y''(t)*dt + 1/2*y'''(t)*dt^2
% yi + 1/2*dt*(y'(t+dt) + y'(t))) + E = y(t) + y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y'''(t)*dt^3 
% 0  + 1/2*dt*(y'(t+dt) + 0)      + E = 0    + 1/2*y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y'''(t)*dt^3 
% 0  + 1/2*dt*(y'(t) + y''(t)*dt + 1/2*y'''(t)*dt^2 + 0)      + E = 0    + 1/2*y'(t)*dt + 1/2*y''(t)*dt^2 + 1/6 y'''(t)*dt^3 
% 0  + 1/2*dt*( + 1/2*y'''(t)*dt^2 + 0)      + E = 0    +  + 1/6 y'''(t)*dt^3
% E = (1/6-1/4) y'''(t)*dt^3 = (2/12-3/12) y'''(t)dt^3 = -1/2 y'''*dt^3

%syms x; d=4; n=4; xx=x.^(0:n), v=vander_1d(xx,n), vd=vanderd_1d(0,n,d), simplify(vd*v^-1)

syms y y1 y2 y3 y4 y5 dt e erel0
y5 = 0
y4 = 0
ynext  = y  + y1*dt + 1/2*y2*dt^2 + 1/6*y3*dt^3 + 1/24*y4*dt^4 + 1/120*y5*dt^5
dynext = y1 + y2*dt + 1/2*y3*dt^2 + 1/6*y4*dt^3 + 1/24*y5*dt^4
eq = ynext == y + dt/2*(y1 + dynext) + e
e = solve(eq,e)

dt0 = solve(abs(e) == e0*abs(y),dt)

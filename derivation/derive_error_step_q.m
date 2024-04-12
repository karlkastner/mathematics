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

syms q

syms y d1y0 d2y0 d3y0 d4y0 d5y0 dt erel0
ynext  = y  + d1y0*dt + 1/2*d2y0*dt^2 + 1/6*d3y0*dt^3 + 1/24*d4y0*dt^4 + 1/120*d5y0*dt^5
d1ynext = d1y0 + d2y0*dt + 1/2*d3y0*dt^2 + 1/6*d4y0*dt^3 + 1/24*d5y0*dt^4
eq = ynext == y + dt*((1-q)*d1y0 + q*d1ynext) + e
e = solve(eq,e)

syms c2 c3;
dt0 = simplify(solve(erel0*abs(y) == abs(c2*dt*dt + c3*dt*dt*dt),dt,'maxdegree',3))
simplify(dt0(1),'ignoreanalyticconstraints',true)

%b = (27*c3^2*e - 2*c2^3 + 3^(1/2)*c3*(4*c2^3*e - 27*c3^2*e^2)^(1/2)*3i)^(1/3)

%- c2/(3*c3) + (2^(2/3)*a)/(6*c3) + (2^(1/3)*c2^2)/(3*c3*a)
%a = (27*c3^2*e - 2*c2^3 + 3^(1/2)*c3*(4*c2^3*e - 27*c3^2*e^2)^(1/2)*3i)^(1/3)
%a = (-27*c3^2*e + 2*c2^3 - 3^(1/2)*c3*(4*c2^3*e - 27*c3^2*e^2)^(1/2)*3i)^(1/3)
if (0)
a = cbrt(27*c3^2*e - 2*c2^3 + 3*c3*sqrt(-12*c2^3*e + 81*c3^2*e^2));
dt0=(cbrtr(2)*(a/2 + c2^2/a) - c2)/(3*c3)
end

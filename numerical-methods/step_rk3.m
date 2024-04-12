% 2023-06-29 11:30:13.466441896 +0200
% Karl Kastner, Berlin
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
%  explicit
%  kutta's third order rungek-kutta method 
%
function [y,e] = step_rk3(t,dt,y0,dy_dt_fun)
	dy0   = dy_dt_fun(t,y0);
	y_mid  = y0 + 0.5*dt*dy0;
	dy_mid = dy_dt_fun(t+0.5*dt,y_mid);
	y1   = y0 - dt*dy0 + (2*dt)*dy_mid;
	dy1  = dy_dt_fun(t+0.5*dt,y_mid);
	y    = y0 + dt*(dy0 + 4*dy_mid + dy1)/6;
	if (nargout()>1)
		y1 = y0 + dt*dy_mid;
		e  = y - y1;
	end
end


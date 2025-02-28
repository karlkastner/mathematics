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
% second order for the reaction equation
%
function [y,e] = step_react_heun(t,dt,y0,dy_dt_fun)
	% TODO for stochastic steps, the random numbers at the second step have to be the same
	dy0 = dy_dt_fun(t,y0);
	% predictor step
	yp  = y0 + dt*dy0;
	%y_  = y.*exp(dt*a./y);
	% corrector-step
	dyp = dy_dt_fun(t+dt,yp);
	y   = y0 + 0.5*dt*(dy0 + dyp);
	if (nargout()>1)
		e = 0.5*dt*(dy0 - dyp);
	end
end


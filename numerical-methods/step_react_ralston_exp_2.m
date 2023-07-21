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
function y = step_react_ralston_exp_2(t,dt,y,dy_dt_fun)
	a   = dy_dt_fun(t,y);
	% predictor step
	%y_  = y + dt*a;
	y_  = y.*exp(dt*a./y);
	% corrector-step
	a_  = dy_dt_fun(t+dt,y_);
	%y   = y + 0.5*dt*(a + a_);
	y   = y.*exp(0.5*dt*(a./y + a_./y_));
end


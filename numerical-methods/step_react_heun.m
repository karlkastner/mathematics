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
%% second order accurate for the reaction equation
%% note: time step has to be adapted externally
function [y,emax,dt_next] = step_react_heun(t,dt,y0,dy_dt_fun,abstol,reltol)
	% TODO for stochastic steps, the random numbers at the second step have to be the same
	dy0 = dy_dt_fun(t,y0);
	% predictor step
	yp  = y0 + dt*dy0;
	%y_  = y.*exp(dt*a./y);
	% corrector-step
	dyp = dy_dt_fun(t+dt,yp);
	y   = y0 + 0.5*dt*(dy0 + dyp);
	if (nargout()>1)
		% The time step is optimal for Euler's method,
		% thus larger than what would be required for Heun's method 
		%
		% y_euler = y0 + dt*dy0;
		% e_euler = 2*abs(y_euler - y_heun)
		%         = 2*dt*abs(dy0 - 1/2*(dy0+dyp))
		%	  =   dt*abs(dy0 - dyp)
		% C dt^2  = e_euler
		% C dt_next^2 = tol
		% dt_next^2/dt^2 = e_euler/tol
		% dt_next = dt*sqrt(tol/e)
		rmsy0 = rms(y0);
		tol   = max(abstol,reltol*rmsy0); 
		e     = 0.5*dt*(dy0 - dyp);
		% note that only the Linf (max) norm is suitable for PDEs
		% as changes might be concentrated in specific locations
		emax    = max(abs(e));
		dt_next = dt*sqrt(tol/emax);
	end
end


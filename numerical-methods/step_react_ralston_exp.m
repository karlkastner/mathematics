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
% second order positivity preserving scheme for the reaction equation
%
% differently to the standard ralston method, the exp is not linearized
% to preserve positivity
function y = step_react_ralston_exp(t,dt,y,afun,bfun)
	% inhomogeneous half-step : prediction
	b  = bfun(t,y);
	y_ = y + 0.5*dt*b;
	% inhomogeneous half-step : correction
	b_ = bfun(t+0.5*dt,y_);
	y  = y + 0.25*dt*(b+b_);

	% homogeneous part : prediction
	a   = afun(t,y);
	y_  = y.*exp(dt*a);
	% corrector-step homogeneous part
	a_  = afun(t+dt,y_);
	y   = y.*exp(0.5*dt*(a + a_));

	% half-step inhomogeneous part
	b  = bfun(t+dt/2,y);
	y_ = y + 0.5*dt*b;
	b_ = bfun(t+dt,y_);
	y = y + 0.25*dt*(b+b_);
end


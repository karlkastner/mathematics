% 2021-06-23 16:03:48.389683349 +0200
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
%
% dy/dt = -z + c_y*(1-h)*y + e_y
% dz/dt = -y + c_z*(1-h)*z + e_z
%
% d2y_dt2 = -dy_dt + c_z*( -dh/dt*z - h*dz_dt) + ez
%
function [t,y] = nonlinear_oscillator_noisy(T,y0,a,c,ca,s,d,dt)
	%ydot = @(t,x) [0,-1; 1, 0]*x + (1-hypot(x(1),x(2))).*x
	%ydot = @(t,x) [0,-1; 1, 0]*x + (1-hypot(x(1),x(2))).*x
	%ydot = @(t,x) [0,-1; 1, 0]*x + c*(1-hypot(x(1),x(2)).^2).*x + s*randn(2,1);
	hat = 1;
	A = [0,-1;1,0];
	% c*(x - x^3 - y^2*x)
        % + c*(x-x.^2) <- simpler, but unstable
	if (0) % == s)
		[t,y] = ode23s(ydot,T,y0);
	else
		%[t,y] = euler(ydot,T,y0);
		[t,y] = ivp_euler_forward2(@ydot,T,cvec(y0),dt);
	end
% only rotation is not possible, bc it depends on phase if e changes amplitude or phase
% if y == +/- 1, dy/dt = z == 0 -> adding to y changes amplitude
% if y == 0, dy/dt = max,  -> adding to y changes phase
function ydot = ydot(t,x)
	e = 2*pi*rand();
	ydot = ( a*[0,-1; 1, 0]*x ... % oscillator
                        +  s*sqrt(dt)*randn(2,1).*[1;0] ... % perturbation
                        ... + s*sqrt(dt)*randn(2,1).*(1 - x.^2) ... % perturbation of mostly phase, cannot be restored without global reference
                        ... + s*sqrt(dt)*randn(2,1).*(x.^2) ... % perturbation of mostly amplitude
			... + 0.01*s*([cos(e)-1,sin(e);-sin(e),cos(e)-1])*x ...
        		+  c*(hat-hypot(x(1),x(2))).*x ... % general recovery (non-linear)
        		+ ca*(hat-hypot(x(1),x(2))).*x.^2.*sign(x) ... % amplitude recovery (non-linear)
        		... + c*(hat-hypot(x(1),x(2))).*x ... % amplitude recovery (non-linear)
			);
	ydot(2) = ydot(2) - d*x(2);
	

end
end


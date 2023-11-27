% 2023-07-17 10:22:21.008600936 +0200
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
function z = step_react_advection_diffusion(t,y,afun,bfun)
	% y' = a*y + b
	% y(t) = (exp(a*t)*(b/a + y0) - b/a);
	% fixed point iteration
	% this is better solved with NR, but this requires solution of a linear system
	kmax = 100;
	kdx = 0;
	o = ones(prod(obj.n),1);
	% pade iteration
	% TODO use function (bcfgf)
	z0 = z;
	a0 = afun(t,z0);
	b0 = bfun(t,z0);
	while (1)
		kdx = kdx+1;
		a = afun(t,z);
		b = bfun(t,z);
		%y = exp(0.5*dt*(a0+a)).*y0 + 0.5*dt*(b + b0);
		z_ = ((1 + 0.5*dt*a0).*y0 + 0.5*dt*(b + b0))./(1 - 0.5*a*dt);
		dz = rms(z - z_);
		z = z_;
		if (dz<tol)
			break;
		end
		if (~isfinite(dz))
			error(['isnan', num2str(kdx), ' ', num2str(dz)]);
		end
		if (kdx>kmax)
			error(['no convergence ', num2str(kdx), ' ', num2str(dz)]);
		end
	end
end


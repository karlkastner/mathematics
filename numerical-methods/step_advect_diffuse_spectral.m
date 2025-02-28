% Mon  2 May 14:18:38 CEST 2022
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
%% solution to the heat equation via Fourier transform
%% with periodic boundary conditions
%%
%% dz/dt = (ex*Dx^2 + ey*Dy^2)*z, z(0) = 0, from 0 to t
%%
%% with fundamental solution:
%% 	z = 1/(4*pi*t)^(n/2)*exp(-r.^2/(4*t))
%%
%% input:
%% t, scalar, : time step
%% dx = [dx,dy] : spatial step size
%% n  = [nx,ny] : number of grid points
%% z0 : nx x ny : initial value
%% e : 2x1 diffusion in each dimension
%% 
%% output:
%%
%% zt : [nx,ny] = value at time t
%%
%% n.B. when dt*e is small and the input is high frequent (discontinuous),
%% then the output is oscillating and not positivity preserving,
%% however, iterative application results in a smooth solution which is exact,
%% i.e. diffuse(100*dt) == 100 x diffuse(dt)
%%
%% function zt = diffuse_spectral(t,z0,dx,e)
%function [zt,fg] = step_advect_diffuse_spectral(dt,dx,a,e,z0,isreal_)
function [zt,fg] = step_advect_diffuse_spectral(z0,dt,a,e,n,L,isreal_)
	ndim = length(L); 
	if (nargin()<7)
		isreal_ = false;
	end	

	% convolve with Gaussian in fourier space
	switch (ndim)
	case {1}
		z0 = cvec(z0);
		%n = length(z0);

		%L = n.*dx;
		ox = 2*pi*fourier_axis(L,n);
		ox = cvec(ox);
		fz0 = fft(z0);
		fg = exp(((-dt*e(1))*ox - (1i*dt*a(1))).*ox);
		fz = fg.*fz0;
		zt = ifft(fz);
	case {2}
		isvector_ = isvector(z0);
		if (isvector_)
			z0 = reshape(z0,n);
		end
		%n = size(z0);
		%L = n.*dx;
		[fx, fy] = fourier_axis_2d(L,n);
		ox = 2*pi*fx;
		oy = 2*pi*fy';
		fz = fft2(z0);
		fg = exp(   ((-dt*e(1))*ox + 1i*dt*a(1)).*ox ...
			  + ((-dt*e(2))*oy + 1i*dt*a(2)).*oy ...
			);
		%fg = exp(-eor2*t -1i*a(1)*ox*dt -1i*a(2)*oy*dt);
		%fg = exp(-(or.*or).*e*t);
		fz = fg.*fz;
		zt = ifft2(fz);
		if (isvector_)
		zt = flat(zt);
		end
		%s = sum(zt(:));
		%zt = zt./s;

	otherwise
		error('not yet implemented');
	end
	if (isreal_)
		zt = real(zt);
	end
end % step_diffuse_spectral


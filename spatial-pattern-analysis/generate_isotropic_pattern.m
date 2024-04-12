% 2021-06-23 21:35:41.688025320 +0200
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% function [z, x, y, Lx, Ly, xx, yy, xe, ye] = generate_isotropic_pattern(fc,n,L,angle0_rad,p,q,scale,st,rotarg,scalearg)
%%
%% spot pattern of unit amplitude
%% output : z : pattern
%%	   x : x-coordinate
%%	   y : y-coordinate
%% 
%% note : rotation, scaling and displacement cannot be fully independently controlled
% s_r,fc_r,order_r
% s_s,rho_s,order_s
% Note : z_gap = 1 - z_spot

% TODO make stoch parameters independent of L and n
% TODO make correlation parameters for rot and scale similar
function [z, x, y, Lx, Ly, xx, yy, xe, ye] = generate_isotropic_pattern(fc,n,L,angle0_rad,p,q,scale,st,rotarg,scalearg)
	if (nargin()<4 || isempty(angle0_rad))
		% rotation of pattern
		angle0_rad=0;
	end
	if (nargin()<5 || isempty(p))
		p = 1;
	end
	if (nargin<6)
		q = 1;
	end
	if (nargin<7)
		scale = true;
	end
	if (nargin()<8)
		st = [];
	end
	if (nargin()<9)
		rotarg = [];
	end
	if (nargin()<10)
		scalearg = [];
	end

	if (scale)
		Lx = L(1);
		s = sqrt(3)/2;
		if (length(L)<2)
			Ly = Lx/s;
		else
			Ly = L(2);
		end
		nx = round(n(1));
		if (length(n)<2)
			ny = round(n(1)/s);
		else
			ny=n(2);
		end
	else
		Lx = L(1);
		Ly = Lx;
		s = 1;
		nx = n(1);
		if (length(n)>1)
			ny = n(2);
		else
			ny = n;
		end
	end

	x = Lx*(double(innerspace(0,1,nx)'-0.5));
	y = Ly*(double(innerspace(0,1,ny)')-0.5);

	xx = repmat(x,1,ny);
	yy = repmat(y,1,nx)';
	%d={0,0,0};
	% stochastic translation (origin)
	if (~isempty(st))
		%xe = brownian_noise_2d_fft(Lx,[nx,ny]);
		%ye = brownian_noise_2d_fft(Ly,[nx,ny]);
		%xx = xx+sbm*xe;
		%yy = yy+sbm*ye;
		nxy = max(nx,ny);
		%if (1)
		ex_t  = brownian_field(0.5,nxy);
		ey_t  = brownian_field(0.5,nxy);
		ex_t  = st*ex_t(1:nx,1:ny); 
		ey_t  = st*ey_t(1:nx,1:ny); 
		%xx = xx+sbm*ex;
		%yy = yy+sbm*ey;
		%d={0,0,0};
		%else	
		%	d = {};
		%	for idx=1:3
		%		d_  = brownian_field(0.5,nxy);
		%		d_  = d_(1:nx,1:ny); 
		%		d{idx} = sbm*d_;
		%	end
		%end
	else
		ex_t = 0;
		ey_t = 0;
	end


	% stochastic rotation (direction)
	if (~isempty(rotarg));
		[ex_r,ey_r,ea] = random_field_rotated(n,L,rotarg{:});
	else
		ex_r = 0;
		ey_r = 0;
	end

	% stochastic scaling (wavelength)
	if (~isempty(scalearg))
		rs = 1;
		[ex_s,ey_s] = random_field_scaled(L,n,scalearg{:});
		%s_s,rs,rho_s,order_s);
	else
		ex_s = 0;
		ey_s = 0;
	end

	% perturb coordinates
	xx = xx + ex_t + ex_r + ex_s;
	yy = yy + ey_t + ey_r + ey_s;

	o = 2*pi*fc;
	z = 0;
	siz = [nx,ny];
	for idx=1:3
		% n.b. when a is randomly perturbed, then there are brighter and darker stripes
		a  = (idx-1)*pi/3+angle0_rad;
		R  = full(rot2(a));
		xy = R*[xx(:),yy(:)]';
		% projected/rotated coordinate
		xp = reshape(xy(1,:),siz);

		% cosine shifted and scaled to [0,1]
		% nb a large perturbation of the amplitude is going to introduce stripes of connected points
		%zi = 1/2*(1+cos(o*x));
		zi = cos(o*xp);
		%(xr+d{idx})));

		% transformation (harmonics)
		%zi = (1-zi.^p).^q;

		% shift and scale back to [-1,1]
		%zi = 2*zi - 1;

		z=z+zi;
	end
	% shift minimum to 0
	z = (2*z+3)/9;

	% transformation (harmonics)
	z = ((1-(1-z).^q).^p);
end % generate_isotropic_stochastic_pattern


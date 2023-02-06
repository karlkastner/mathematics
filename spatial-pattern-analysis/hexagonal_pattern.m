% 2021-06-23 21:35:41.688025320 +0200
%
% function [z, x, y, xx, yy, xe, ye] = hexagonal_pattern(fc,n,L,a0,scale,sbm,p,q)
%
% spot pattern of unit amplitude
% output : z : pattern
%	   x : x-coordinate
%	   y : y-coordinate
% 
% Note : z_gap = 1 - z_spot
function [z, x, y, Lx, Ly, xx, yy, xe, ye] = hexagonal_pattern(fc,n,L,a0,scale,sbm,p,q)
	if (nargin()<4 || isempty(a0))
		% rotation of pattern
		a0=0;
	end
	if (nargin()<5 || isempty(scale))
		scale = true;
	end
	if (nargin()<6)
		sbm = [];
	end
	if (nargin<7)
		p = 1;
	end
	if (nargin<8)
		q = 1;
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
	if (1)
	if (~isempty(sbm))
		xe = brownian_noise_2d_fft(Lx,[nx,ny]);
		ye = brownian_noise_2d_fft(Ly,[nx,ny]);
		xx = xx+sbm*xe;
		yy = yy+sbm*ye;
	end
	end

	o = 2*pi*fc;
	z= 0;
	siz = [nx,ny];
	for idx=1:3
		% n.b. when a is randomly perturbed, then there are brighter and darker stripes
		a  = (idx-1)*pi/3+a0;
		R  = full(rot2(a));
		xy = R*[xx(:),yy(:)]';
		% n.b when xr is randomly perturbed then there are also brighter and darker stripes
		xr  = reshape(xy(1,:),siz);
 		%yr  = reshape(xy(2,:),siz);

		if (0)
		if (~isempty(sbm))
			xe = brownian_noise_2d_fft(L,[nx,ny]);
			xr = xr+sbm*xe;
			if (idx==2)
			ye = xe;
			end
		else
			xe = 1;
			ye = 0;
		end
		end

		% cosine shifted and scaled to [0,1]
		% nb a large perturbation of the amplitude is going to introduce stripes of connected points
		zi = 1/2*(1+cos(o*xr));

		% transformation (harmonics)
		zi = (1-zi.^p).^q;

		% shift and scale back to [-1,1]
		zi = 2*zi - 1;

		z=z+zi;
	end
	% shift minimum to 0
	z = (3-2*z)/9;

end % hexagonal pattern


% 2021-06-23 21:35:41.688025320 +0200
% spot pattern of unit amplitude
% output : z
% 
% Note : z_gap = 1 - z_spot
function z = hexagonal_pattern(f,n,L,a0,scale)
	p = 1;
	q = 1;
	if (nargin()<4)
		a0=0;
	end
	if (nargin()<5)
		scale = true;
	end

	if (scale)
		s = sqrt(3)/2;
		Lx = L;
		Ly = L/s;
		nx = round(n(1));
		ny = round(n(1)/s);
	else
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
	o = 2*pi*f;
	xx = repmat(x,1,ny);
	yy = repmat(y,1,nx)';
	z= 0;

	for idx=1:3
		a  = (idx-1)*pi/3+a0;
		R  = full(rot2(a));
		xy = R*[xx(:),yy(:)]';
		x  = reshape(xy(1,:),size(x0));
 		y  = reshape(xy(2,:),size(y0));

		% cosine shifted and scaled to [0,1]
		zi = 1/2*(1+cos(o*x));

		% transformation (harmonics)
		zi = (1-zi.^p).^q;

		% shift and scale back to [-1,1]
		zi = 2*zi - 1;

		z=z+zi;
	end
	% shift minimum to 0
	z = (3-2*z)/9;

end % hexagonal pattern


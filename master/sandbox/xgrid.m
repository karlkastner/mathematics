% Wed Nov  2 18:52:01 MSK 2011
% Karl Kästner, Berlin

% computes an irregular 1D-grid with more points located at the given
% centrepoint depending on the choosen order
%
% for two output arguments grid return the unscaled grid and a scale factor
% this can be used to generate integer-based difference matrices without round
% off error and subsequent apllication of the scale factor
%
% n : numer of grid points
% L : half domain size, domain : [-L ... L]
% p : grid spacing power, can be odd
% p = 1 : equal spacing (regular grid), h_const = 1/n
% p > 1 : gridpoints closer to the centrepoint with h_min ~ (1/n)^p
function [X scale] = xgrid(n,p,x0,L,vargrid)
	if (nargin < 2)
		p = 1.0;
	end
	if (nargin < 3)
		x0 = 0.0;
	end
	if (x0/L > 1 || x0/L < -1)
		'error'
		return
	end
	if (nargin < 4)
		L = 1.0;
	end
	if (nargin < 5)
		vargrid = 'quadratic';
	end
	n = n-2;
	switch (vargrid)
		case {'quadratic'}
			x0 = x0/L;
			h  = ((1-x0)^(1/p) + (1+x0)^(1/p))/(n+1);
			xl = (x0-(1+x0)^(1/p))/h;
			xr = (x0+(1-x0)^(1/p))/h;
			X  = [xl:1:xr];
			d  = X-x0/h;
			X  = x0/h^p + sign(d).*abs(d.^p);
			scale = L*h^p;
		% integer central
			%X = (1:n/2+1).^2;
			%X = [-fliplr(X) X]
			%scale = L/X(end);
			X_bc = X';
		case {'inverse'}
			X = (0:n+1)/(n+1) - 0.5;
			Y = 0.5./X;
			X = -L*fftshift(Y)/(n+1);
			X_bc = X';
		case {'logarithmic'}
			h=1/(n+1);
			c = 1;
			b = (h)^(2/n);
			X = b.^(0:n/2);
			X = fliplr(X);
			X = L*[-fliplr(X) X]/X(end);
			X_bc = X';
		case {'central'}
			X = L*[0:n/4 (1:n/2)/(n/2+1) n/4+1:n/2+1]'/(n/2+1) - L/2;
			h = 1/(n+1);
			X_bc = L*h*[0:2:n/2 h*((n*(n+1)/2+1):2:((n/2+1)*(n+1)-1)) n/2+1:2:n+1]' - L/2;
		case {2,'2'} 
			% second order accurate grid
			X = sum1(1:n/2+1) - 0.5;
			X_bc = L*[-fliplr(X) X]/X(end);
			X_bc = X_bc';
		case {4,'4'}
			% fourth order accurate grid
			X = sum3(1:n/2+1) - 0.5;
			X_bc = L*[-fliplr(X) X]/X(end);
			X_bc = X_bc';
		otherwise 
			% do not use x = [ 0.5*h 1-0.5h ], as bc is only valid for x = [h 1-h]
			% first order accurate grid
			X_bc = 2*(0:n+1)' - (n+1);
			scale = L/(n+1);
	end
	% scale
	if (nargout < 2)
		X = scale*X;
	end
end % function grid


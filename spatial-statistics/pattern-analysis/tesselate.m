% 2025-02-18 17:00:23.199551587 +0100 tile.m
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
function [z,x,dx,dy] = tile(z1,L1,sdx,sL,s)
	n1 = size(z1,1);
	if (nargin()<2)
		L1 = n1;
	end
	if (nargin()<3)
		sdx = 1;
	end
	if (nargin()<2)
		sL = 1;
	end

	%n1  = size(z1,1);
	%L1  = n1(1)/2;
	%dx1 = L1/n1(1);

	n   = sdx*sL*n1; 
	L   = sL*L1;

	x1 = ((0:n1-1)'+1/2)*L1/n1;
	x  = ((0:n-1)'+1/2)*L/n;
	y = x';	

if (0)	
	dx = brownian_field(0.5,n);
	dy = brownian_field(0.5,n);

	[fx,fy,fr] = fourier_axis_2d([L,L],[n,n]);
	fdx = fft2(dx);
	fdy = fft2(dy);
	fdx(fr>0.5) = 0;
	fdy(fr>0.5) = 0;
	dx = ifft2(fdx);
	dy = ifft2(fdy);
	
	sx = 3;
	sy = 3;
else
	% low frequent noise
	[fx,fy,fr] = fourier_axis_2d([L,L],[n,n]);
	%S = lowpass2d_continuous_pdf(fr,0.1,1);
	S = lowpass2d_discrete_pdf([L,L],[n,n],1,s);
	S = S/(sum(S,'all')*(fx(2)-fx(1))^2);
	dx = randn(n);
	dy = randn(n);
	dx = ifft2(S.*fft2(dx));
	dy = ifft2(S.*fft2(dy));
	sx = s;
	sy = s;
end

	xx = x + sx*dx;
	yy = y + sy*dy;

	% padd 3x3
	z1 = [z1,z1,z1;
  	    z1,z1,z1;
	      z1,z1,z1];
	x1 = [x1-L1;x1;x1+L1];
	z = interp2(x1,x1,z1,mod(xx,L1),mod(yy,L1),'linear');


% Tue 28 Feb 17:25:12 CET 2023
% Karl Kastner, Berlin
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
% this is not grid averaged
function [a,x,y,a0,z0] = geometric_ar1_2d_generate(lmu,lsd,theta,L,n)
	x = ((0:n-1)'/n-0.5)*L;
	x = ifftshift(x);
	y = x;
	%acfun = @(x1,x2,y1,y2) exp(-theta*hypot(x2-x1,y2-y1));
	%ac = acfun(0,x,0,y');
	acfun = @(x,y) exp(-theta*hypot(x,y));
	ac = acfun(x,y');
	%ac = ifftshift(ac);
	% generate correlated normals
	z0 = randn(n);
	S  = real(fft2(ac));
	df = 1/L;
	S  = S/(sum(S,'all')*df^2);
	T  = sqrt(S);
	z  = ifft2(T.*fft2(z0));
	z  = z/std(z,[],'all');
	a  = exp(lmu+lsd*z);
	z0 = z0/std(z0,[],'all');
	a0 = exp(lmu+lsd*z0);
end


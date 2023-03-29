% Mon 13 Mar 14:54:28 CET 2023
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
%
function [R,x,y] = anisotropic_pattern_acf(L,n,f0,sxy)
	x = fourier_axis(n(1)/L(1),n(1));
	y = fourier_axis(n(2)/L(2),n(2));
	k0 = 2*pi*f0;
	R = cos(k0*rvec(x)).*exp(-pi*k0*hypot(sxy(1)^2*rvec(x),sxy(2)^2*cvec(y)));
end


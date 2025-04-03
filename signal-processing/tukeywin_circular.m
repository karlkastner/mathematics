% 2025-02-28 11:08:34.054188412 +0100 mathematics/signal-processing/tukeywin_circular.m
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
%% circular Tukey window in two dimensions
%
% TODO make elliptic for rectangular domains
% TODO aliasing (integer weights at boundary)
function w= tukeywin_circular(n,pw)
	[fx,fy,frr] = fourier_axis_2d([1,1],n);
	f2 = n(1)/2;
	f1 = (1-pw)*f2;
	w = (1+cos(pi*(frr - f1)/(f2 - f1)))/2;
	%w = (1+cos(pi*(frr - pw*f2)/(f2 - pw*f2)))/2;
	w(frr < f1) = 1;
	w(frr > f2) = 0; % > or >=?
end


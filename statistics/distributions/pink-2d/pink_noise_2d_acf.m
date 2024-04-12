% Fri  5 Jan 10:31:23 CET 2024
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% note that this seems to be independent of dx?
% TODO consider L
function [R,x,y,rr] = pink_noise_2d_acf(n,L)
	[x,y,rr] = real_axis_2d(n,L);
	%y = (0:n(2)-1)-n(1)/2;
	%r = hypot(x,y);

	R = 0;
	%N = length(r);
	for k=1:n(1)
		R  = R + besselj(0,2*pi*k*rr/n(1))./k;
	end
	R = R/R(1);
end


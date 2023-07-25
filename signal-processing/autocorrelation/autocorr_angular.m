% Sat 14 Jan 15:02:43 CET 2023
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
function [Rt,ti] = autocorr_angular(R,L,m)
	n = size(R);
	[x,y,r,t] = fourier_axis_2d(L,n);
	Lmax = min(L)/2;
	ti   = innerspace(-pi,pi,m);
	dt   = ti(2)-ti(1);
	id   = floor((t + pi)/dt)+1;
	fdx  = abs(r)<Lmax;
	Rt   = accumarray(id(fdx),R(fdx),[length(ti),1],@mean, 0);
end


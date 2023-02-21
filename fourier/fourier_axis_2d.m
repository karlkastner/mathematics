% Wed 11 May 13:10:56 CEST 2022
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
%% function [fx, fy, fr, ft, Tx, Ty, mask, N] = fourier_axis_2d(L,n)
function [fx, fy, fr, ft, Tx, Ty, mask, N] = fourier_axis_2d(L,n)
	[fx,Tx]=fourier_axis(L(1),n(1));
	[fy,Ty]=fourier_axis(L(2),n(2));
	fr = hypot(cvec(fx),rvec(fy));
	ft = atan2(rvec(fy),cvec(fx));
end


% Wed 27 Apr 11:02:07 CEST 2022
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
%% truncated, not wrapped at the end
function [R,x,y] = lowpass2d_discrete_acf(L,n,La)
	if (length(La) == 1)
		La(2) = La(1);
	end
	[x,y]    = fourier_axis_2d(n./L,n);
	r_div_a  = hypot(cvec(x/La(1)),rvec(y/La(2)));
	% autocorrelation
	R = exp(-r_div_a);
end



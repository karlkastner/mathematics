% 2018-08-29 14:27:21.089395412 +0200
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
%% integrate y along x with the trapezoidal rule
function cumint_y = int_trapezoidal(x,y)
	if (isvector(y))
		y = cvec(y);
	end
	if (isvector(x))
		x = cvec(x);
	end
	dx = diff(x);
	cumint_y = [zeros(1,size(y,2)); 0.5*cumsum((y(1:end-1,:)+y(2:end,:)).*dx)];
end


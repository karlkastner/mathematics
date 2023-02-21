% 2023-01-15 17:54:09.272681131 +0100
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
function [fc,Sc] = gamma_mode(a,b,flag)
	if (a<1)
		fc = 0;
	else
	%	fc = a-1/b;
		fc = (a-1)*b;
	end
	if (nargin()<3 || ~flag)
		Sc = gampdf(fc,a,b);
	else
		Sc = gampdf_man(fc,a,b)
	end
end


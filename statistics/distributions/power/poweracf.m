% 2025-03-05 11:42:01.992567208 +0100
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

% S = s*1/f^a
% R = s*2 int_0^infty 1/f^a exp(-2 i pi f x) df 
%   ~     1/x^(1-a)

function acf = poweracf(x,a,dim)
	if (nargin()<3)
		dim = 1;	
	end
	scale = (2*pi)^(a-1);
	acf   = scale./(x.^(1+a));
end


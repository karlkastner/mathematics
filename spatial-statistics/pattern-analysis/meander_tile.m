% 2025-02-18 14:47:09.386638589 +0100
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
function [z1] = meander_tile(n1)

	% prototype tile
	z1 = zeros(n1,n1);
	dir = [ 0,1
                1,0
                0,-1
               -1,0];
                
	z1(1,1:n1-1) = 1;
	id = [1,n1-1];
	ddx = 1;
	n1_ = n1(1)-2;
	while (n1_>0)
		ddx = mod(ddx,4)+1;
		id_(1) = id(1)+n1_*dir(ddx,1);
		id_(2) = id(2)+n1_*dir(ddx,2);
		if (dir(ddx,1)~=0)
			z1(id(1):dir(ddx,1):id_(1),id(2)) = 1;
		else
			z1(id(1),id(2):dir(ddx,2):id_(2)) = 1;
		end
		id = id_;
		n1_ = n1_ - 1;
	end
end


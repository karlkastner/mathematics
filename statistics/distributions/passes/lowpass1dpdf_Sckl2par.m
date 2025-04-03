% 2025-03-29 10:38:49.558579511 +0100
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
function [a,b] = lowpass2d_Sckl2par(Sc,kl,par0)
	if (nargin<3)
	par0 = [1,1];
	end
	par = lsqnonlin(@(p) lowpass1dpdf([0,kl],p(1),p(2)) - Sc*[1,0.5],par0);
	a = par(1);
	b = par(2);
end


% 2023-03-03 21:02:02.088631189 +0100
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
% only valid for fr<1
function Sb = lowpass2d_continuous_pdf(fr,a,order,m)
	n = m;
	Sb = 0;
	or = 2*pi*fr;
	for l=0:n
		Sb = Sb + (-1)^l*(1/2)^(2*l)/factorial(l).^2*or.^(2*l)*a^-(2*n+2)*factorial(2*l+1);
	end
	Sb = 2*pi*Sb;
end


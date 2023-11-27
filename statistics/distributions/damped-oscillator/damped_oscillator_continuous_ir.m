% Fri 27 Oct 10:17:52 CEST 2023
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
function [ir] = damped_oscillator_continuous_ir(x,a1,a2)
	r = abs(a1)/(2*a2);
	fc = sqrt(4*a2 - a1^2)./(2*a2);
	ir = (a1*x>0).*exp(-abs(x)*r).*sin(fc*x)./sqrt(4*a2 - a1^2);
end


% Fri 27 Oct 09:32:14 CEST 2023
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
function [fc,Sc]=damped_oscillator_continuous_mode(a1,a2)
	fc=  sqrt(a2 - a1.^2/2)./(2*pi*a2);
	Sc = 4*a1*(4*a2.^2)./(4*a2.*a1.^2 - a1^4);
end


% Fri 27 Oct 09:59:21 CEST 2023
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
function [a1,a2] = damped_oscillator_continuous_mode2par(fc0,Sc0)
	poly = [64*Sc0^2*fc0^6*pi^6, + 16*Sc0^2*fc0^4*pi^4, - 4*Sc0^2*fc0^2*pi^2 + 32, - Sc0^2];
	a2 = roots(poly);
	% determine real root
	[~,fdx] = min(abs(imag(a2)));
	a2 = real(a2(fdx));
	a1 = sqrt(2*a2*(1-4*a2*fc0^2*pi^2));
end


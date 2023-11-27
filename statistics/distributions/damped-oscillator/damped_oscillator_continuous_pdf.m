% Thu 26 Oct 17:30:42 CEST 2023
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
function S = damped_oscillator_continuous_pdf(f,a1,a2)
	o = 2*pi*f;
	%S = 1./(1 + a1*1i*o - a2*o.^2).*1./(1 - a1*1i*o - a2*o.^2);
        S  = 1./(1 + (a1^2 - 2*a2)*o.^2 + a2^2*o.^4);
	% normalization
	S = 4*a1*S;
end


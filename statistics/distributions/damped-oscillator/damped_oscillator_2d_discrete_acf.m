% Wed 25 Oct 11:04:13 CEST 2023
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
% autocorrelation function of the discrete phase noise integration process
% (damped oscillator with phase noise)
% 	(z - a1 dz/dx + a2 d^2z/dx^2) == e
function [acf1,S1,AA] = damped_oscillator_1d_discrete_acf(n,L,a1,a2,a3)
	[S1,AA] = damped_oscillator_2d_discrete_pdf(n,L,a1,a2,a3);
	% note : the result is real, but imaginary part has to be stripped due to round off errors
	acf1 = real(ifft2(S1));
	acf1 = acf1/acf1(1,1);
end


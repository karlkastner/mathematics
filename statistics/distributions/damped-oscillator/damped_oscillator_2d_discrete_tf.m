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
% transfer function of the discrete phase noise integration process
% (damped oscillator with phase noise)
% 	(z - a1 dz/dx + a2 d^2z/dx^2) == e
function [tf,A] = damped_oscillator_2d_discrete_tf(n,L,a1,a2,a3)
	[ir,A] = damped_oscillator_2d_discrete_ir(n,L,a1,a2,a3);
	tf = fft2(ir);
end


% Fri 10 Nov 11:14:47 CET 2023
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
% t : vector with angles between 0 and 2 pi-dt
% S : spectral density (strictly positive)
% piso : 1 if isotropic, 0 if ansisotropic
function [piso,t0,c,St_] = isisotropic(t,St)
	n = length(t);
	t = cvec(t);
	A = [ones(n,1),sin(2*t),cos(2*t)];
	c = A \ cvec(St);
	caniso = hypot(c(2),c(3));
	ciso   = c(1)-caniso;
	piso   = (0.5*caniso+ciso)./(caniso + ciso);
	t0 = atan2(c(3),c(2));
	St_ = A*c;
end


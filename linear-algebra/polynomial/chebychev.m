% Fri  7 Apr 14:02:38 CEST 2017
% Karl Kastner, Berlin
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
%% c.f. Dronkers 1964, eq. 8.15, p. 300
%% chebycheff polynomials
%% function c = chebychev(x,n)
function c = chebychev(x,n)
	if (issym(x))
		syms c
	else
		c=zeros(1,n);
	end
	c(1) = 1;
	if (n>0)
		c(2) = x;
	end
	for idx=2:n
		c(idx+1) = 2*x*c(idx) - c(idx-1);
	end
end


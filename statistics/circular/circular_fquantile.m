% Wed 28 Jun 11:29:06 CEST 2023
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
function x0 = circular_fquantile(fun,q,nx)
	if (nargin()<3)
		nx = 100;
	end
	x  = linspace(-pi,pi,nx)';
	p  = fun(x);
	c  = cumsum(p);
	c  = c/c(end);
	% make unique
	fdx = [true;c(2:end)~=c(1:end-1)];

	x0 = interp1(c(fdx),x(fdx),q,'linear');
end


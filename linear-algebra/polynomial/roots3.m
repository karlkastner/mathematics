% 2020-03-27 03:14:57.400149719 +0100
% Karl Kästner, Berlin
%
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
%
function r = roots3(c)
	d0 = c(:,2).^2 - 3*c(:,1).*c(:,3);
	d1 = 2*c(:,2).^3 - 9*c(:,1).*c(:,2).*c(:,3) + 27*c(:,1).^2.*c(:,4);
	% TODO sign choice inside sqrt so that C is not zero
	C  = cbrt(0.5*(d1 + sqrt(d1.^2 - 4*d0.^3)));
	e = (0.5*(-1 + sqrt(-3))).^(0:2);
	r = -1./(3*c(:,1)).*(c(:,2) + C*e + d0./(C*e));
end


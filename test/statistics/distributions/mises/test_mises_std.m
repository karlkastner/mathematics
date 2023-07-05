% 2023-06-28 11:12:42.929899604 +0200

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

n = 1e4;
x = innerspace(-pi,pi,n+1);
k = 1;
y  = misespdf(x,0,k);
dx = x(2)-x(1);
s2 = sum(y.*x.^2*dx)
sd = sqrt(s2)

s2_ = mises_var(0,k)
sd_ = mises_std(0,k)
sd/sd_



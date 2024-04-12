% 2024-01-04 19:52:41.137057243 +0100
% Karl Kastner, Berlin
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
function xd = downsample1(x)
	xd          = 0.5*x(1:2:end-1,:);
	xd(1,:)     = xd(1,:) + 0.25*(x(end,:) + x(2,:));
	xd(2:end,:) = xd(2:end,:) + 0.25*(x(2:2:end-1,:) + x(4:2:end,:));
end



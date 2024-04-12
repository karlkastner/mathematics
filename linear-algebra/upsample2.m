% 2024-01-04 19:42:58.456358721 +0100
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
%
function xr = upsample2(x)
	n             = 2*size(x,2);
	% for allocation, element n is assigned first
	xr(:,n)       = 0.5*(x(:,1)+x(:,end));
	xr(:,1:2:n-1) = x;
	xr(:,2:2:n-2) = 0.5*(x(:,1:end-1)+x(:,2:end));
end


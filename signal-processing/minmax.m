% Tue 13 Dec 15:43:42 CET 2022
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
function [mindx,maxdx] = minmax(x)
	maxdx_=1;
	mindx_=1;
	mindx = 1;
	maxdx = 1;
	% find min
	for idx=1:length(x)
		if (x(idx) < x(mindx_))
			mindx_ = idx;
			maxdx_ = idx;
		elseif (x(idx) > x(maxdx_))
			maxdx_ = idx;
			if ((x(maxdx_)-x(mindx_)) > (x(maxdx)-x(mindx)))
				maxdx = maxdx_;
				mindx = mindx_;
			end
		end
	end
end


% 2023-06-22 11:41:13.438050211 +0200
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% note: conventionally, S(0) is weighted by 1/2
function S = periodogram_normalize(S,df,fullaxis)
	if (fullaxis)
		% full axis
		S = S/(sum(S)*df);
	else
		% half axis
		S = 2*S/(sum(S)*df);
	end
end


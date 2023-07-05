% 2023-06-22 11:40:39.233529966 +0200
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
% note S(x,0) and S(0,y) are conventionally weighted by 1/2
%      S(0,0) by 1/4
function S =  periodogram_normalize_2d(S,df)
	S = 2*S./(sum(S,'all')*df(1)*df(2));
end


% Tue  5 Dec 12:22:23 CET 2023
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
% hellinger distance between two probability distributions
%
% hd = 1/2 int w (S1^(1/2) - S2^(1/2))^2 dx
%      1/2 (int w S1 dx - 2 int w sqrt(S1*S2) + int w S2 dx)
%    = 1 -  int w sqrt(S1 cdot S2)^2 dx
%
% 0 <= hd <= 1
% 
function hd = hellinger_distance(S1,S2,df,w)
	S1 = cvec(S1);
	S2 = cvec(S2);
	if (nargin()>3)
		w = cvec(w);
		% normalize
		S1 = S1./(sum(w.*S1)*df);
		S2 = S2./(sum(w.*S2)*df);
		hd = 1 - sum(w.*sqrt(S1.*S2))*df; 
	else % w = 1
		hd = 1 - sum(sqrt(S1.*S2))*df; 
	end
end


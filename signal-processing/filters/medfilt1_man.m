% Di 26. Mai 08:49:39 CEST 2015
% Karl Kastner, Berlin
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
%% moving median filter, supports columnwise operation
function [Y, S, L, U] = medfilt1_man(X,n)
	Y = zeros(size(X));
	S = zeros(size(X));
	L = zeros(size(X));
	U = zeros(size(X));
	l = floor(n/2);
	r = round(n/2)-1;
	for idx=1:size(X,1)
		[Y(idx,:), S(idx,:), L(idx,:), U(idx,:)] = median_man( ...
			X(max(1,idx-l):min(end,idx+r),:)');
	end
end % medfilt1_man


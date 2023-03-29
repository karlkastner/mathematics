% 2021-06-23 20:22:48.970336838 +0200
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
%% this function is computationally inefficient and serves merely for illustration
%% and tests
function x = lowpass2d_2(x,L)
	s = 5;
	if (L>0)
	xf = -s*L:s*L;
	f = exp(-1/2*(xf/L).^2)
	f = f/sum(f);
	x = conv2(x,f,'same');
	x = conv2(x,f','same');
	%x = x-conv2(x,f,'same');
	end
end

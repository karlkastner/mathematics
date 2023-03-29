% 2023-03-13 14:34:20.327295481 +0100
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
function [S,x,y,R] = anisotropic_pattern_pdf(L,n,f0,sxy)
	[R,x,y] = anisotropic_pattern_acf(L,n,f0,sxy);
	S = real(ifft2(R));
	df = 1./L;
	S = 2*S./(sum(S,'all')*df(1)*df(2));
end


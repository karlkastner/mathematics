% 2021-10-21 09:29:25.938016430 +0200
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
%
%% integrate the spectral density over the positive half axis
%
function int_S = spectral_density_area(fx,S)
	if (isvector(S))
		S = cvec(S);
	end
	fdx = (fx>=0);
	int_S = int_trapezoidal(fx(fdx),S(fdx,:));
end


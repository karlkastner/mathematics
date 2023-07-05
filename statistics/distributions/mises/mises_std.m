% Wed 28 Jun 11:03:59 CEST 2023
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
% circular standard deviation
function cs1 = mises_cstd(mu,k)
	cs2 = mises_var(mu,k);
	R  = 1 - cs2;
	cs1 = sqrt(-2*log(R));
end


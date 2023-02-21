% Sa 1. Aug 13:16:03 CEST 2015
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
%% convert kendall rank correlation coefficient to the person product moment
%% correlation coefficient
%%
%% c.f. Kruskal, 1958, p. 823
function rho = kendall_to_pearson(tau)
	rho = sin(pi*tau/2);
end


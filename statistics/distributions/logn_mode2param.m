% 2022-03-30 12:46:31.571547062 +0200
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution

function [lmu,lsd] = logn_mode2param(xm,ym)
	% sqrt(2*pi)*lc*Sc = exp(-1/2*s^2)/s
	lsd2 = lambertw(1./(2*pi*xm.^2.*ym.^2));
	lsd = sqrt(lsd2);
	lmu = log(xm)+lsd2;
end


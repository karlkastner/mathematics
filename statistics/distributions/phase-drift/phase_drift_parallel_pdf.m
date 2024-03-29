% Sat 11 Jun 15:06:15 CEST 2022
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
% function S = spectral_density_brownian_phase_across(fy,sy)
%
function S = spectral_density_brownian_phase_across(fy,sy)
	if (issym(fy) || issym(sy))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	ky = 2*pi_*fy;
	S = 4*sy.^2./(sy.^4+ky.*ky);
end

